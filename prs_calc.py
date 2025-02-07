import sys
import pandas as pd
from PyQt6.QtWidgets import QApplication, QWidget, QVBoxLayout, QPushButton, QLabel, QFileDialog, QTextEdit


class PRSApp(QWidget):
    def __init__(self):
        super().__init__()
        self.gwas_data = {}
        self.snp_data = []
        self.initUI()
        self.load_gwas_data()

    def initUI(self):
        self.setWindowTitle("Polygenic Risk Score Calculator")
        self.setGeometry(200, 200, 600, 400)
        layout = QVBoxLayout()

        self.snp_label = QLabel("SNP List: Not loaded yet")
        layout.addWidget(self.snp_label)

        self.load_snp_button = QPushButton("Load SNP List")
        self.load_snp_button.clicked.connect(self.load_snp_file)
        layout.addWidget(self.load_snp_button)

        self.calculate_button = QPushButton("Calculate PRS")
        self.calculate_button.clicked.connect(self.calculate_prs)
        layout.addWidget(self.calculate_button)

        self.result_text = QTextEdit()
        self.result_text.setReadOnly(True)
        layout.addWidget(self.result_text)

        self.setLayout(layout)

    def load_gwas_data(self):
        try:
            gwas_file = "gwas_data.tsv"  # GWAS file stored within the program
            gwas_df = pd.read_csv(gwas_file, sep='\t', low_memory=False, dtype=str)
            gwas_df.dropna(subset=['SNPS', 'OR or BETA'], inplace=True)
            gwas_df['OR or BETA'] = pd.to_numeric(gwas_df['OR or BETA'], errors='coerce')
            gwas_df.dropna(subset=['OR or BETA'], inplace=True)
            self.gwas_data = dict(zip(gwas_df['SNPS'], gwas_df['OR or BETA']))
        except Exception as e:
            print(f"Error loading GWAS data: {e}")

    def load_snp_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select SNP List", "", "Excel Files (*.xlsx *.xls)")
        if file_path:
            self.snp_label.setText(f"SNP List: {file_path}")
            try:
                df = pd.read_excel(file_path)
                self.snp_data = df.dropna(subset=['rsID', 'Genotype']).to_dict(orient='records')
                print("SNP data loaded successfully.")
            except Exception as e:
                print(f"Error reading SNP file: {e}")
                self.snp_label.setText("Error loading SNP file.")

    def calculate_prs(self):
        if self.snp_data:
            prs_score = 0
            genotype_mapping = {'AA': 0, 'AB': 1, 'BB': 2}  # Example mapping, adjust as needed
            for snp in self.snp_data:
                rsid = snp.get('rsID')
                genotype = snp.get('Genotype')
                if rsid in self.gwas_data and genotype in genotype_mapping:
                    try:
                        prs_score += genotype_mapping[genotype] * self.gwas_data[rsid]
                    except ValueError as e:
                        print(f"Error in genotype value for {rsid}: {e}")
                else:
                    print(f"Missing data for rsID: {rsid} or genotype: {genotype}")
            interpretation = self.interpret_prs(prs_score)
            result_text = f"PRS Score: {prs_score:.4f}\n\n{interpretation}"
            self.result_text.setText(result_text)
        else:
            self.result_text.setText("Please load your SNP file.")

    def interpret_prs(self, score):
        if score < -1:
            return "You are in a genetically low-risk group. However, living a better life is always a good option."
        elif -1 <= score < 1:
            return "You are in an average genetic risk group. A healthy lifestyle is recommended."
        else:
            return "You are in a genetically high-risk group. Consult your doctor for precautions."


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = PRSApp()
    window.show()
    sys.exit(app.exec())

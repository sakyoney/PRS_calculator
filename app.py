from __future__ import annotations

import uuid
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd
from flask import Flask, flash, redirect, render_template, request, send_file, url_for
from werkzeug.utils import secure_filename

BASE_DIR = Path(__file__).resolve().parent
UPLOAD_DIR = BASE_DIR / "uploads"
RESULT_DIR = BASE_DIR / "results"
ALLOWED_EXTENSIONS = {"csv", "tsv", "txt", "xls", "xlsx"}

UPLOAD_DIR.mkdir(exist_ok=True)
RESULT_DIR.mkdir(exist_ok=True)

app = Flask(__name__)
app.config["SECRET_KEY"] = "prs-demo-secret-key"
app.config["MAX_CONTENT_LENGTH"] = 16 * 1024 * 1024  # 16 MB


def allowed_file(filename: str) -> bool:
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    if suffix == ".csv":
        return pd.read_csv(path)
    return pd.read_csv(path, sep=None, engine="python")


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(c).strip() for c in out.columns]
    return out


def pick_column(df: pd.DataFrame, candidates: Iterable[str]) -> Optional[str]:
    norm_map = {str(col).strip().lower(): col for col in df.columns}
    for candidate in candidates:
        if candidate.lower() in norm_map:
            return norm_map[candidate.lower()]
    return None


def genotype_to_dosage(genotype: str, effect_allele: str) -> Optional[float]:
    if pd.isna(genotype):
        return None
    g = str(genotype).strip().replace(" ", "")
    if not g:
        return None

    for sep in ["/", "|", "\\", "-"]:
        if sep in g:
            parts = [p.upper() for p in g.split(sep) if p]
            if len(parts) == 2:
                return float(sum(1 for allele in parts if allele == effect_allele.upper()))

    if len(g) == 2 and g.isalpha():
        parts = [g[0].upper(), g[1].upper()]
        return float(sum(1 for allele in parts if allele == effect_allele.upper()))

    try:
        return float(g)
    except ValueError:
        return None


def compute_prs(gwas_df: pd.DataFrame, sample_df: pd.DataFrame) -> tuple[pd.DataFrame, float, int]:
    gwas_df = normalize_columns(gwas_df)
    sample_df = normalize_columns(sample_df)

    gwas_cols = {
        "snp": pick_column(gwas_df, ["snp", "rsid", "markername", "variant", "id"]),
        "beta": pick_column(gwas_df, ["beta", "effect", "weight", "log_or", "estimate"]),
        "effect_allele": pick_column(gwas_df, ["effect_allele", "ea", "a1", "allele1", "risk_allele"]),
    }
    sample_cols = {
        "snp": pick_column(sample_df, ["snp", "rsid", "markername", "variant", "id"]),
        "dosage": pick_column(sample_df, ["dosage", "effect_allele_dosage", "ea_dosage", "score"]),
        "genotype": pick_column(sample_df, ["genotype", "gt", "geno"]),
    }

    missing_gwas = [key for key, value in gwas_cols.items() if value is None]
    if missing_gwas:
        raise ValueError(f"GWAS file is missing required columns: {', '.join(missing_gwas)}")

    if sample_cols["snp"] is None:
        raise ValueError("Sample file must include a SNP/rsid column.")
    if sample_cols["dosage"] is None and sample_cols["genotype"] is None:
        raise ValueError("Sample file must include either a dosage column or a genotype column.")

    gwas_work = gwas_df[[gwas_cols["snp"], gwas_cols["beta"], gwas_cols["effect_allele"]]].copy()
    gwas_work.columns = ["SNP", "BETA", "EFFECT_ALLELE"]
    gwas_work["SNP"] = gwas_work["SNP"].astype(str).str.strip()
    gwas_work["EFFECT_ALLELE"] = gwas_work["EFFECT_ALLELE"].astype(str).str.strip().str.upper()
    gwas_work["BETA"] = pd.to_numeric(gwas_work["BETA"], errors="coerce")
    gwas_work = gwas_work.dropna(subset=["SNP", "BETA", "EFFECT_ALLELE"])

    sample_keep = [sample_cols["snp"]]
    if sample_cols["dosage"]:
        sample_keep.append(sample_cols["dosage"])
    if sample_cols["genotype"]:
        sample_keep.append(sample_cols["genotype"])
    sample_work = sample_df[sample_keep].copy()
    rename_map = {sample_cols["snp"]: "SNP"}
    if sample_cols["dosage"]:
        rename_map[sample_cols["dosage"]] = "DOSAGE"
    if sample_cols["genotype"]:
        rename_map[sample_cols["genotype"]] = "GENOTYPE"
    sample_work = sample_work.rename(columns=rename_map)
    sample_work["SNP"] = sample_work["SNP"].astype(str).str.strip()

    merged = pd.merge(gwas_work, sample_work, on="SNP", how="inner")
    if merged.empty:
        raise ValueError("No overlapping SNPs were found between the GWAS and sample files.")

    if "DOSAGE" in merged.columns:
        merged["DOSAGE"] = pd.to_numeric(merged["DOSAGE"], errors="coerce")

    if "GENOTYPE" in merged.columns:
        genotype_dosage = merged.apply(
            lambda row: genotype_to_dosage(row.get("GENOTYPE"), row["EFFECT_ALLELE"]), axis=1
        )
        if "DOSAGE" not in merged.columns:
            merged["DOSAGE"] = genotype_dosage
        else:
            merged["DOSAGE"] = merged["DOSAGE"].fillna(genotype_dosage)

    merged = merged.dropna(subset=["DOSAGE", "BETA"])
    if merged.empty:
        raise ValueError("Overlapping SNPs were found, but none had usable dosage/genotype values.")

    merged["PRS_COMPONENT"] = merged["DOSAGE"] * merged["BETA"]
    prs_score = float(merged["PRS_COMPONENT"].sum())
    snp_count = int(merged.shape[0])

    result_cols = ["SNP", "EFFECT_ALLELE", "BETA", "DOSAGE", "PRS_COMPONENT"]
    if "GENOTYPE" in merged.columns:
        result_cols.insert(4, "GENOTYPE")
    result_df = merged[result_cols].sort_values("SNP")
    return result_df, prs_score, snp_count


@app.route("/", methods=["GET"])
def index():
    return render_template("index.html")


@app.route("/calculate", methods=["POST"])
def calculate():
    gwas_file = request.files.get("gwas_file")
    sample_file = request.files.get("sample_file")

    if not gwas_file or not sample_file:
        flash("Please upload both a GWAS file and a sample file.")
        return redirect(url_for("index"))

    if not allowed_file(gwas_file.filename) or not allowed_file(sample_file.filename):
        flash("Unsupported file type. Allowed: csv, tsv, txt, xls, xlsx.")
        return redirect(url_for("index"))

    request_id = uuid.uuid4().hex[:10]
    gwas_path = UPLOAD_DIR / f"{request_id}_gwas_{secure_filename(gwas_file.filename)}"
    sample_path = UPLOAD_DIR / f"{request_id}_sample_{secure_filename(sample_file.filename)}"
    result_path = RESULT_DIR / f"prs_result_{request_id}.csv"

    try:
        gwas_file.save(gwas_path)
        sample_file.save(sample_path)

        gwas_df = read_table(gwas_path)
        sample_df = read_table(sample_path)
        result_df, prs_score, snp_count = compute_prs(gwas_df, sample_df)
        result_df.to_csv(result_path, index=False)

        preview = result_df.head(20).to_dict(orient="records")
        return render_template(
            "index.html",
            prs_score=round(prs_score, 6),
            snp_count=snp_count,
            preview=preview,
            download_name=result_path.name,
        )
    except Exception as exc:  # noqa: BLE001
        flash(f"Calculation failed: {exc}")
        return redirect(url_for("index"))


@app.route("/download/<filename>", methods=["GET"])
def download(filename: str):
    file_path = RESULT_DIR / filename
    if not file_path.exists():
        flash("Result file not found.")
        return redirect(url_for("index"))
    return send_file(file_path, as_attachment=True)


if __name__ == "__main__":
    app.run(debug=True)

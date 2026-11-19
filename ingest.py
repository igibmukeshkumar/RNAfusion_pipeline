import os
import re
import csv
import json
import argparse
from pathlib import Path
from lxml import etree
from openai import OpenAI

# ============================================================
# Config
# ============================================================
FUSION_TERMS = re.compile(
    r"\b("
    r"gene fusion|fusion gene|fusion transcript|fusion transcripts|"
    r"fusion RNA|fusion mRNA|chimeric transcript|chimeric transcripts|"
    r"chimeric RNA|chimeric mRNA|fusion protein|fusion proteins|"
    r"chimeric protein|chimeric proteins|gene rearrangement|"
    r"gene rearrangements|genomic rearrangement|genomic rearrangements|"
    r"chromosomal rearrangement|chromosomal rearrangements|"
    r"translocation|translocations|breakpoint|breakpoints|"
    r"transcript joining|exon.{0,30}junction|exon.{0,30}fusion|"
    r"[A-Z][A-Z0-9]{1,15}::[A-Z][A-Z0-9]{1,15}|"
    r"[A-Z][A-Z0-9]{1,15}--[A-Z][A-Z0-9]{1,15}|"
    r"[A-Z][A-Z0-9]{1,15}-[A-Z][A-Z0-9]{1,15}|"
    r"[A-Z][A-Z0-9]{1,15}/[A-Z][A-Z0-9]{1,15}"
    r")\b", re.I
)

# ============================================================
# NEMOTRON PROMPT
# ============================================================
SYSTEM_PROMPT = """
You are an expert biomedical literature curator.

Here your job is curate reliable fusion curation.

RECALL: Extract EVERY distinct supported fusion pair. If X is fused with Z, Y, A, B and C, 
return all five pairs. Never stop at the first.

Reject genes merely mentioned together, gene lists, expression/co-occurrence,
rearrangement without fusion evidence, unsupported claims, and engineered
fusion proteins unrelated to biological gene fusion.


Return ONLY valid JSON:
{
  "gene_a": "",
  "gene_b": "",
  "fusion_type": "GENE_FUSION|RNA_FUSION|FUSION_TRANSCRIPT|FUSION_PROTEIN|REARRANGEMENT_ONLY|NOT_FUSION|UNCERTAIN",
  "confirmed": true,
  "confidence": 0.0,
  "reliability": "HIGH|MEDIUM|LOW",
  "evidence": "",
  "reason": ""
}
"""

# ============================================================
# XML extraction
# ============================================================

def text(node):
    return " ".join("".join(node.itertext()).split())

def pubmed_id(article, name):
    x = article.xpath(f'.//*[local-name()="{name}"][1]')
    return text(x[0]) if x else ""

def pubmed_articles(root):
    return root.xpath('.//*[local-name()="PubmedArticle"]')

# article have article and name both, as PMID is also supplied here 
def read_pubmed(article):
    pmid = pubmed_id(article, "PMID")
    pmc = ""
    ids = article.xpath('.//*[local-name()="ArticleId"]')
    for x in ids:
        typ = (x.get("IdType") or "").lower()
        if typ in {"pmc", "pmcid"}:
            pmc = text(x)
            break

    title = article.xpath('.//*[local-name()="ArticleTitle"][1]')
    abstract = article.xpath('.//*[local-name()="AbstractText"]')
    chunks = ([text(title[0])] if title else []) + [text(x) for x in abstract]
    return pmid, pmc, [x for x in chunks if x]

def read_pmc(root):
    article = root.xpath('.//*[local-name()="article"][1]')[0]
    pmid = ""
    pmc = ""

    for x in article.xpath('.//*[local-name()="article-id"]'):
        typ = (x.get("pub-id-type") or "").lower()
        value = text(x)
        if typ == "pmid":
            pmid = value
        elif typ in {"pmcid", "pmc"}:
            pmc = value

    title = article.xpath('.//*[local-name()="article-title"][1]')
    abstract = article.xpath('.//*[local-name()="abstract"]')
    paragraphs = article.xpath('.//*[local-name()="body"]//*[local-name()="p"]')

    chunks = []
    if title:
        chunks.append(text(title[0]))
    chunks += [text(x) for x in abstract]
    chunks += [text(x) for x in paragraphs]

    return pmid, pmc, [x for x in chunks if x]

def read_xml(path):
    root = etree.parse(str(path)).getroot()

    if pubmed_articles(root):
        for article in pubmed_articles(root):
            yield read_pubmed(article)
    else:
        yield read_pmc(root)

# ============================================================
# Fusion retrieval
# ============================================================

def fusion_passages(chunks):
    hits = []
    for i, chunk in enumerate(chunks):
        if not FUSION_TERMS.search(chunk):
            continue
        passage = " ".join(chunks[max(0, i - 1):min(len(chunks), i + 2)])
        if passage not in hits:
            hits.append(passage)
    return hits

# ============================================================
# Nemotron
# ============================================================

def ask_nemotron(client, model, passage):
    response = client.chat.completions.create(
        model=model,
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": f"PASSAGE:\n{passage}"},
        ],
        temperature=0.0,
        max_tokens=500,
        response_format={"type": "json_object"},
        extra_body={"chat_template_kwargs": {"enable_thinking": False}},
    )
    result = response.choices[0].message.content.strip()
    result = re.sub(r"^```json\s*|\s*```$", "", result)
    try:
        return json.loads(result)
    except Exception:
        return {}

def parse_nemotron(x):
    gene_a = str(
        x.get("gene_a") or x.get("gene1") or x.get("gene_1") or ""
    ).strip().upper()
    gene_b = str(
        x.get("gene_b") or x.get("gene2") or x.get("gene_2") or ""
    ).strip().upper()

    if not gene_a or not gene_b or gene_a == gene_b:
        return None

    return {
        "gene_a": gene_a,
        "gene_b": gene_b,
        "fusion": f"{gene_a}--{gene_b}",
        "fusion_type": str(x.get("fusion_type", "")),
        "confidence": float(x.get("confidence", 0.0) or 0.0),
        "confirmed": bool(x.get("confirmed", False)),
        "reliability": str(x.get("reliability", "")),
        "evidence": str(x.get("evidence", "")),
    }

# ============================================================
# HGNC
# ============================================================

def load_hgnc(path):
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames or []
        col = next(
            (
                x for x in fields
                if x.lower().replace(" ", "_")
                in {"symbol", "approved_symbol", "gene_symbol", "hgnc_symbol"}
            ),
            fields[0],
        )
        return {
            row[col].strip().upper()
            for row in reader
            if row.get(col)
        }

def hgnc_status(gene_a, gene_b, hgnc):
    a = gene_a.upper() in hgnc
    b = gene_b.upper() in hgnc
    if a and b:
        return "BOTH_MATCHED", "YES"
    if a or b:
        return "ONE_MATCHED", "NO"
    return "NONE_MATCHED", "NO"

# ============================================================
# Main
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description="Nemotron-only RNA/gene fusion literature curation"
    )
    parser.add_argument("-i", action="append", required=True)
    parser.add_argument("--gene-map", required=True)
    parser.add_argument("--nvidia-threshold", type=float, default=0.5)
    parser.add_argument(
        "--nvidia-model",
        default="nvidia/nemotron-3.5-lightning-30b-a3b"
    )
    parser.add_argument(
        "--nvidia-base-url",
        default="https://integrate.api.nvidia.com/v1"
    )
    parser.add_argument(
        "--skip-processed",
        action="store_true",
        help="Skip article IDs already present in the output TSV"
    )
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    client = OpenAI(
        api_key=os.environ.get("NVIDIA_API_KEY", "local"),
        base_url=args.nvidia_base_url,
    )
    hgnc = load_hgnc(args.gene_map)

    files = []
    for directory in args.i:
        files.extend(sorted(Path(directory).glob("*.xml")))

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "id", "pmid", "pmc", "gene_a", "gene_b", "fusion",
        "fusion_type", "nvidia_confidence", "validated_with_hgnc",
        "hgnc_approved", "nvidia_confirmed", "reliability",
        "evidence", "source",
    ]

    processed = set()
    if args.skip_processed and out.exists():
        with open(out, encoding="utf-8", newline="") as f:
            processed = {
                row["id"]
                for row in csv.DictReader(f, delimiter="\t")
                if row.get("id")
            }

    write_header = not out.exists() or out.stat().st_size == 0
    mode = "a" if args.skip_processed else "w"

    article_count = 0
    skipped_count = 0

    with open(out, mode, encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        if write_header or mode == "w":
            writer.writeheader()

        for n, path in enumerate(files, 1):
            try:
                for pmid, pmc, chunks in read_xml(path):
                    article_count += 1
                    article_id = pmc or pmid or path.stem

                    if args.skip_processed and article_id in processed:
                        skipped_count += 1
                        continue

                    article_rows = []

                    for passage in fusion_passages(chunks):
                        x = parse_nemotron(
                            ask_nemotron(client, args.nvidia_model, passage)
                        )
                        if x is None:
                            continue

                        confirmed = (
                            x["confirmed"]
                            and x["confidence"] >= args.nvidia_threshold
                        )

                        hgnc_result, approved = hgnc_status(
                            x["gene_a"], x["gene_b"], hgnc
                        )

                        article_rows.append({
                            "id": article_id,
                            "pmid": pmid,
                            "pmc": pmc,
                            "gene_a": x["gene_a"],
                            "gene_b": x["gene_b"],
                            "fusion": x["fusion"],
                            "fusion_type": x["fusion_type"],
                            "nvidia_confidence": round(x["confidence"], 3),
                            "validated_with_hgnc": hgnc_result,
                            "hgnc_approved": approved,
                            "nvidia_confirmed": "YES" if confirmed else "NO",
                            "reliability": x["reliability"],
                            "evidence": x["evidence"],
                            "source": str(path),
                        })

                    for row in article_rows:
                        writer.writerow(row)
                    f.flush()

                    # Only mark articles with output rows as processed.
                    if article_rows:
                        processed.add(article_id)

            except Exception as e:
                print(f"[ERROR] {path}: {e}")

            if n % 100 == 0:
                print(
                    f"Files: {n}/{len(files)} | Articles: {article_count} | "
                    f"Skipped: {skipped_count}"
                )

    print(f"Done: {out}")
    print(f"Articles scanned: {article_count}")
    print(f"Articles skipped: {skipped_count}")


if __name__ == "__main__":
    main()

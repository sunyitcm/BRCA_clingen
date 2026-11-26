"""Return cell values from the exact locked Table 9 workbook for human review."""
import re
from pathlib import PurePosixPath
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from .core import R
from .net import getj
from .reference import file_sha

MAIN = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def attach_table_notes(report, directory):
    rows = [c["source"] for c in report["candidates"] if c["source"]["table"] == "T9"]
    rows.extend(r for r in report["facts"].get("workbook_rows_without_assignment", [])
                if r["table"] == "T9")
    if not rows:
        return
    lock = next(x for x in getj(R / "rules/attachments.lock.json")
                if x["gene"] == report["gene"] and x["name"].startswith("Specifications_Table9"))
    path = directory / lock["name"]
    if file_sha(path) != lock["sha256"]:
        raise ValueError("Table 9 工作簿与冻结版本不符")
    by_row = {r["row"]: r for r in rows}
    found = set()
    with ZipFile(path) as book:
        wb = ET.fromstring(book.read("xl/workbook.xml"))
        sheet = next(s for s in wb.iter(MAIN+"sheet") if s.get("name") == rows[0]["sheet"])
        relationships = ET.fromstring(book.read("xl/_rels/workbook.xml.rels"))
        target = next(r.get("Target") for r in relationships if r.get("Id") == sheet.get(REL+"id"))
        member = target.lstrip("/") if target.startswith("/") else str(PurePosixPath("xl") / target)
        strings = []
        if "xl/sharedStrings.xml" in book.namelist():
            strings = ["".join(t.text or "" for t in si.iter(MAIN+"t"))
                       for si in ET.fromstring(book.read("xl/sharedStrings.xml"))]
        with book.open(member) as handle:
            for _, node in ET.iterparse(handle, events=("end",)):
                if node.tag != MAIN+"row":
                    continue
                row_id = int(node.get("r"))
                if row_id in by_row:
                    cells = {}
                    for cell in node.findall(MAIN+"c"):
                        column = re.sub(r"[0-9]+$", "", cell.get("r"))
                        value = cell.findtext(MAIN+"v", "")
                        if cell.get("t") == "s":
                            value = strings[int(value)]
                        elif cell.get("t") == "inlineStr":
                            value = "".join(t.text or "" for t in cell.iter(MAIN+"t"))
                        cells[column] = value
                    row = by_row[row_id]
                    if (cells.get("A"), cells.get("B")) != (report["gene"], row["hgvsc"]):
                        raise ValueError("Table 9 物理行身份不符")
                    # Match the actual assignment columns; citation text is not an assignment.
                    if (cells.get("D"),cells.get("E"))==(row.get("code"),row.get("strength")) and row.get("code") in ("PS3","BS3") and row.get("strength")=="Strong":
                        row["checked_assignment"]={"gene":cells["A"],"hgvsc":cells["B"],"code":cells["D"],"strength":cells["E"],"cells":f"D{row_id}:E{row_id}","workbook_sha256":lock["sha256"]}
                    row["review_values"] = {c: cells.get(c) for c in "FGHIJKLMN"}
                    found.add(row_id)
                node.clear()
                if found == by_row.keys():
                    break
    if found != by_row.keys():
        raise ValueError("Table 9 物理行缺失")

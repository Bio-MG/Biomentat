#!/usr/bin/env python3
"""tools/status_ledger.py — Canonical status ledger generator and C6-C10 guards.

Subcommands / modes:
  --build : Reads docs/ROADMAP.md, folder indexes, and scan-inventory profile;
            writes canonical JSON registry to docs/agent/status_ledger.json
            and updates the <!-- LEDGER:BEGIN --> ... <!-- LEDGER:END --> block
            in docs/ROADMAP.md deterministically (idempotent).
  --check : Validates docs/agent/status_ledger.json and docs/ROADMAP.md against
            guards C6 (closed enum), C7 (TERMINATED iff archive/), C8 (mandatory
            anchor & deliverable/partial proof), C9 (population > 0), and
            C10 (ROADMAP.md <-> status_ledger.json agreement).

Exit codes:
  0 = All C6-C10 guards pass (or --build succeeded on valid tree/fixture)
  1 = One or more C6-C10 guard violations detected (--check)
  2 = Instrument failure (missing arguments/files or unexpected error)
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True

import argparse
import json
import os
import re
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SKILL_SCRIPTS = (
    _REPO_ROOT / "docs" / "agent" / "skills" / "scan-inventory" / "scripts"
)
if str(_SKILL_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SKILL_SCRIPTS))

import census as _census  # noqa: E402

CORE_FAMILIES = ("proposals", "ideation", "contracts", "archive", "benchmarks")

PREFIX_MAP = {
    "proposals": "PROP",
    "ideation": "IDEA",
    "benchmarks": "BENCH",
    "contracts": "CONTRACT",
    "decisions": "ADR",
    "archive": "ARCH",
    "(root)": "ROOT",
}

KIND_MAP = {
    "proposals": "proposal",
    "ideation": "ideation",
    "benchmarks": "benchmark",
    "contracts": "contract",
    "decisions": "decision",
    "archive": "archive",
    "(root)": "roadmap",
}


def strip_ledger_block(card_text: str, begin_m: str, end_m: str) -> str:
    pattern = re.compile(
        re.escape(begin_m) + r".*?" + re.escape(end_m) + r"\n?",
        re.S,
    )
    return pattern.sub("", card_text)


def make_item_id(family: str, rel: Path) -> str:
    prefix = PREFIX_MAP.get(family, family.upper())
    stem = re.sub(r"[^A-Za-z0-9_-]+", "_", rel.stem).strip("_")
    if stem.upper().startswith(prefix + "_") or stem.upper().startswith(prefix + "-"):
        return stem
    return f"{prefix}-{stem}"


def infer_item_kind(family: str, role: str, rel: Path) -> str:
    if role == "evidence" or family in ("benchmarks", "audits"):
        return "evidence"
    if role == "frozen" or family in ("contracts", "decisions"):
        return "decision"
    upper = rel.stem.upper()
    if "FIX" in upper or "BUG" in upper or "CLEANUP" in upper or "HYGIENE" in upper:
        return "fix"
    if "FUNCTION" in upper:
        return "function"
    if any(k in upper for k in ("FEATURE", "PACKS", "GEMS", "SEURAT", "CCC", "DRAWER")):
        return "feature"
    return "improvement"


def extract_evidence(
    text: str,
    card_row: str = "",
    index_row: str = "",
) -> str:
    fields = {m.group(1): m.group(2).strip() for m in _census.FIELD_RE.finditer(card_row)}
    if fields.get("EVIDENCE"):
        return fields["EVIDENCE"]
    if index_row:
        pm = _census.PROPOSAL_REF_RE.search(index_row)
        if pm:
            return pm.group(0)
        em = _census.EVIDENCE_RE.search(index_row)
        if em:
            return em.group(0)
    if card_row:
        cm = _census.EVIDENCE_RE.search(card_row)
        if cm:
            return cm.group(0)
    tm = _census.EVIDENCE_RE.search(text)
    if tm:
        return tm.group(0)
    return ""


def extract_progress(
    text: str,
    card_row: str,
    state: str,
    notes: list[str],
) -> str:
    fields = {m.group(1): m.group(2).strip() for m in _census.FIELD_RE.finditer(card_row)}
    if fields.get("JALONS"):
        return fields["JALONS"]
    n_done = len(_census.JALON_DONE_RE.findall(text))
    n_open = len(_census.JALON_OPEN_RE.findall(text))
    if n_done + n_open > 0:
        return f"{n_done}/{n_done + n_open}"
    if state == "PARTIAL":
        if any(n.startswith("status_map:") for n in notes):
            return "1/2"
        if "card_hint=PARTIAL" in notes or "keyword+evidence" in notes:
            return "1/2"
    if state == "DELIVERED":
        return "1/1"
    return "0/0"


def extract_next_action(card_row: str, state: str) -> str:
    fields = {m.group(1): m.group(2).strip() for m in _census.FIELD_RE.finditer(card_row)}
    if fields.get("NEXT"):
        return fields["NEXT"]
    if state == "TERMINATED":
        return "none (archived)"
    if state == "FROZEN":
        return "frozen (requires new GO)"
    if state == "DELIVERED":
        return "M3 archive"
    if state == "UNTOUCHED":
        return "arbitrage operateur"
    return "see card"


def resolve_roadmap_anchor(
    fam: str,
    rel: Path,
    registered: bool,
    not_judged: bool,
    card_rel_posix: str,
    card_text: str,
    indexes: dict[str, tuple[str, str]],
) -> str:
    if not registered and not not_judged:
        return ""
    if fam in ("proposals", "benchmarks", "(root)"):
        if _census.whole_token_in(card_text, rel.name):
            # Identify subsection in card if possible
            pos = card_text.find(rel.name)
            if pos != -1:
                prefix_text = card_text[:pos]
                headings = re.findall(r"^###\s+([^\n]+)", prefix_text, re.M)
                if headings:
                    sec = headings[-1].split("—")[0].strip()
                    return f"{card_rel_posix} §1 {sec}"
            return f"{card_rel_posix} §1"
        return "" if not not_judged else f"{card_rel_posix} §1"
    if fam in indexes:
        return indexes[fam][0]
    if fam == "decisions":
        return f"{card_rel_posix} §1 Socle (docs/decisions/)"
    return f"{card_rel_posix} §1"


def scan_and_build_ledger(
    root: Path,
    card: Path,
    raw_card_text: str,
    profile: dict,
) -> tuple[dict, list[str]]:
    markers = profile.get("ledger_markers", {})
    begin_m = markers.get("begin", "<!-- LEDGER:BEGIN -->")
    end_m = markers.get("end", "<!-- LEDGER:END -->")
    card_text = strip_ledger_block(raw_card_text, begin_m, end_m)

    docs_root = root / profile.get("roots", {}).get("docs", "docs")
    families_cfg = profile.get("families", {})
    states = list(profile.get("states", []))
    allowed_enum = set(states)

    if not docs_root.is_dir():
        raise RuntimeError("docs root not found: " + str(docs_root))

    by_family: dict[str, list[Path]] = {fam: [] for fam in families_cfg}
    for d in sorted(docs_root.iterdir()):
        if d.is_dir() and not d.name.startswith("."):
            by_family.setdefault(d.name, [])
    for p in sorted(docs_root.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(root)
        if rel.name.endswith((".bak", ".tmp", ".log", ".pyc")) or rel.name == "status_ledger.json":
            continue
        fam = rel.parts[1] if len(rel.parts) > 2 else "(root)"
        by_family.setdefault(fam, []).append(rel)

    card_paths = set(_census.BACKTICK_PATH_RE.findall(card_text))
    for _cp in list(card_paths):
        card_paths.add(Path(_cp).name)
    indexes: dict[str, tuple[str, str]] = {}
    for fam, cfg in families_cfg.items():
        idx = cfg.get("index")
        if idx:
            indexes[fam] = (idx, _census.safe_read(root / idx))

    try:
        card_rel_posix = card.relative_to(root).as_posix()
    except ValueError:
        card_rel_posix = "docs/ROADMAP.md"

    surface_paths = {c.get("index", "") for c in families_cfg.values() if c.get("index")}
    surface_paths.update({
        "docs/ROADMAP.md",
        "docs/STATUS.md",
        "docs/CONVENTIONS.md",
        "docs/README.md",
        "docs/agent/HANDOFF.md",
        "docs/agent/DECISIONS.md",
        "docs/agent/status_ledger.json",
        card_rel_posix,
    })

    budget_exceeded_paths: set[str] = set()
    for bkey, limit in profile.get("budgets", {}).items():
        if "#" not in bkey:
            bpath = root / bkey
            if bpath.is_file():
                n_lines = len(_census.safe_read(bpath).rstrip().splitlines())
                if n_lines > limit:
                    budget_exceeded_paths.add(bkey)
        elif bkey.endswith("#0"):
            rel_p = bkey.split("#", 1)[0]
            bpath = root / rel_p
            if bpath.is_file():
                txt = _census.safe_read(bpath)
                m = re.search(r"^## 0\.", txt, re.M)
                if m:
                    nxt = re.search(r"^## ", txt[m.end():], re.M)
                    body = txt[m.start(): (m.end() + nxt.start()) if nxt else len(txt)]
                    if len(body.rstrip().splitlines()) > limit:
                        budget_exceeded_paths.add(rel_p)
        elif bkey.endswith("#header"):
            rel_p = bkey.split("#", 1)[0]
            bpath = root / rel_p
            if bpath.is_file():
                txt = _census.safe_read(bpath)
                m = re.search(r"^---\s*$", txt, re.M)
                if m and len(txt[:m.start()].rstrip().splitlines()) > limit:
                    budget_exceeded_paths.add(rel_p)

    def registered_in(fam: str, rel: Path) -> bool:
        if rel.name in card_paths or _census.whole_token_in(card_text, rel.name):
            return True
        if fam not in indexes:
            return False
        body = indexes[fam][1]
        stem = rel.stem
        alts = {rel.name, stem, re.sub(r"_(CONTRACT|PROPOSAL)$", "", stem)}
        return any(_census.whole_token_in(body, alt) for alt in alts)

    items: list[dict] = []
    surfaces: list[dict] = []
    findings: list[str] = []

    for fam, files in sorted(by_family.items()):
        cfg = families_cfg.get(fam, {})
        topology = cfg.get("topology", "grouped")
        fam_allowed = set(cfg.get("allowed_states", states))
        role = cfg.get("role", "work")
        index_text = indexes.get(fam, ("", ""))[1]
        is_not_judged_fam = bool(cfg.get("not_judged") or topology == "grouped")

        for rel in files:
            rel_posix = rel.as_posix()
            text = _census.safe_read(root / rel)
            clean_text = _census.strip_code_fences(text)
            card_row = _census.card_row_for(card_text, rel.name)
            is_b_exceeded = rel_posix in budget_exceeded_paths

            if rel_posix in surface_paths:
                surfaces.append({
                    "path": rel_posix,
                    "family": fam,
                    "flags": ["stale"] if is_b_exceeded else [],
                })
                continue

            nested = len(rel.parts) > 3 and fam != "(root)"

            if is_not_judged_fam and fam not in ("archive", "decisions"):
                tok_m = _census.TOKEN_RE.search(card_row) or _census.INLINE_STATE_RE.search(clean_text)
                if tok_m and tok_m.group(1) == "TERMINATED" and "archive" not in rel.parts:
                    findings.append(f"C7: TERMINATED_OUTSIDE_ARCHIVE: `{rel_posix}`")
                if fam == "(root)" and registered_in(fam, rel):
                    state, notes = _census.classify_state(
                        text, rel, cfg, profile, card_row=card_row, index_text=index_text
                    )
                    fl = _census.flags_for(
                        state, text, True, profile, card_row=card_row, notes=notes, budget_exceeded=is_b_exceeded
                    )
                    ev = extract_evidence(text, card_row)
                    prog = extract_progress(text, card_row, state, notes)
                    nxt = extract_next_action(card_row, state)
                    anchor = resolve_roadmap_anchor(
                        fam, rel, True, True, card_rel_posix, card_text, indexes
                    )
                    items.append({
                        "id": make_item_id(fam, rel),
                        "path": rel_posix,
                        "family": fam,
                        "kind": KIND_MAP.get(fam, fam),
                        "item_kind": infer_item_kind(fam, role, rel),
                        "state": state,
                        "progress": prog,
                        "evidence": ev,
                        "next_action": nxt,
                        "flags": fl,
                        "roadmap_anchor": anchor,
                    })
                continue

            if fam == "archive" and nested:
                continue

            registered = True if fam == "decisions" else registered_in(fam, rel)
            if "archive" not in rel.parts and _census.INLINE_BANNER_RE.search(clean_text):
                findings.append(f"C8: SINGLE_SURFACE_VIOLATION: `{rel_posix}` carries inline state banner/token")

            # Check raw token before classify_state normalizes it
            raw_tok_m = _census.TOKEN_RE.search(card_row) or _census.INLINE_STATE_RE.search(clean_text)
            raw_tok_state = raw_tok_m.group(1) if raw_tok_m else ""

            state, notes = _census.classify_state(
                text, rel, cfg, profile, card_row=card_row, index_text=index_text
            )
            _, idx_row = (
                _census.index_row_for(index_text, rel.name, cfg.get("status_column", "Status"))
                if index_text
                else ("", "")
            )

            # C6 checks
            if any(n.startswith("out_of_enum_state=") for n in notes) or (
                raw_tok_state and raw_tok_state not in allowed_enum
            ):
                bad_st = raw_tok_state or state
                findings.append(f"C6: OUT_OF_ENUM_STATE: `{rel_posix}` state={bad_st}")
            elif state not in allowed_enum or state not in fam_allowed:
                findings.append(f"C6: STATE_NOT_ALLOWED: `{rel_posix}` state={state} family={fam}")

            # C7 checks
            if state == "TERMINATED" and "archive" not in rel.parts:
                findings.append(f"C7: TERMINATED_OUTSIDE_ARCHIVE: `{rel_posix}`")
            if "archive" in rel.parts and state != "TERMINATED":
                findings.append(f"C7: ARCHIVE_NOT_TERMINATED: `{rel_posix}` state={state}")

            # C8 checks
            if "partial_without_closed_jalon->ACTIVE" in notes:
                findings.append(f"C8: PARTIAL_WITHOUT_CLOSED_JALON: `{rel_posix}`")
            if state == "DELIVERED" and not _census.has_deliverable_evidence(text, card_row, idx_row):
                findings.append(f"C8: DELIVERED_WITHOUT_EVIDENCE: `{rel_posix}`")
            if topology in ("direct", "folder_index") and not cfg.get("not_judged") and not registered and not nested:
                owner = indexes[fam][0] if fam in indexes else cfg.get("owner", "card")
                findings.append(f"C8: UNREGISTERED_MISSING_ANCHOR: `{rel_posix}` — owner={owner}")

            fl = _census.flags_for(
                state,
                text,
                registered,
                profile,
                card_row=card_row,
                index_row=idx_row,
                notes=notes,
                budget_exceeded=is_b_exceeded,
            )
            ev = extract_evidence(text, card_row, idx_row)
            prog = extract_progress(text, card_row, state, notes)
            nxt = extract_next_action(card_row, state)
            anchor = resolve_roadmap_anchor(
                fam, rel, registered, bool(cfg.get("not_judged")), card_rel_posix, card_text, indexes
            )

            items.append({
                "id": make_item_id(fam, rel),
                "path": rel_posix,
                "family": fam,
                "kind": KIND_MAP.get(fam, fam),
                "item_kind": infer_item_kind(fam, role, rel),
                "state": state,
                "progress": prog,
                "evidence": ev,
                "next_action": nxt,
                "flags": fl,
                "roadmap_anchor": anchor,
            })

    # C9 population > 0 check
    populations = {f: len(v) for f, v in sorted(by_family.items())}
    empty_fams = [f for f, cnt in populations.items() if cnt == 0]
    for cf in CORE_FAMILIES:
        if populations.get(cf, 0) == 0 and cf not in empty_fams:
            empty_fams.append(cf)
    if empty_fams:
        findings.append("C9: EMPTY_POPULATION: " + ", ".join(sorted(set(empty_fams))))
    if not items:
        findings.append("C9: EMPTY_LEDGER_ITEMS: 0 classified items")

    items.sort(key=lambda x: x["path"])
    surfaces.sort(key=lambda x: x["path"])

    by_state_counts = {st: 0 for st in states}
    by_fam_state: dict[str, dict[str, int]] = {}
    for it in items:
        st = it["state"]
        fm = it["family"]
        if st in by_state_counts:
            by_state_counts[st] += 1
        by_fam_state.setdefault(fm, {s: 0 for s in states})
        if st in by_fam_state[fm]:
            by_fam_state[fm][st] += 1

    ledger = {
        "schema_version": 1,
        "card": card_rel_posix,
        "populations": populations,
        "summary": {
            "total_items": len(items),
            "by_state": by_state_counts,
            "by_family": {k: by_fam_state[k] for k in sorted(by_fam_state)},
        },
        "surfaces": surfaces,
        "items": items,
    }
    return ledger, findings


def render_ledger_block(ledger: dict, profile: dict) -> str:
    markers = profile.get("ledger_markers", {})
    begin_m = markers.get("begin", "<!-- LEDGER:BEGIN -->")
    end_m = markers.get("end", "<!-- LEDGER:END -->")
    states = list(profile.get("states", ["ACTIVE", "PARTIAL", "UNTOUCHED", "FROZEN", "DELIVERED", "TERMINATED"]))

    by_fam = ledger.get("summary", {}).get("by_family", {})
    by_state = ledger.get("summary", {}).get("by_state", {})
    total_items = ledger.get("summary", {}).get("total_items", 0)

    lines = [
        begin_m,
        "### Synthèse générée du registre de statut (`docs/agent/status_ledger.json`)",
        "",
        "> Bloc généré par `python -B tools/status_ledger.py --build` — ne pas éditer à la main (garde `C10`).",
        "",
        "| Famille | Total | " + " | ".join(states) + " |",
        "|---|---:|" + "|".join(["---:"] * len(states)) + "|",
    ]
    for fam in sorted(by_fam.keys()):
        counts = by_fam[fam]
        row_tot = sum(counts.get(s, 0) for s in states)
        cells = [str(counts.get(s, 0)) for s in states]
        lines.append(f"| {fam} | {row_tot} | " + " | ".join(cells) + " |")

    tot_cells = [f"**{by_state.get(s, 0)}**" for s in states]
    lines.append(f"| **Total** | **{total_items}** | " + " | ".join(tot_cells) + " |")

    action_items = [
        it for it in ledger.get("items", [])
        if it.get("family") in ("proposals", "(root)") and it.get("state") in ("ACTIVE", "PARTIAL")
    ]
    if action_items:
        lines += [
            "",
            "#### Chantiers `ACTIVE` / `PARTIAL` sur la carte (`proposals` & `(root)`)",
            "",
            "| ID | État | Avancement | Preuve | Prochaine action | Drapeaux |",
            "|---|---|---|---|---|---|",
        ]
        for it in action_items:
            ev = it.get("evidence") or "—"
            prog = it.get("progress") or "—"
            nxt = it.get("next_action") or "—"
            fl = ",".join(it.get("flags", [])) or "—"
            lines.append(
                f"| `{it['id']}` | {it['state']} | {prog} | {ev} | {nxt} | {fl} |"
            )

    lines.append(end_m)
    return "\n".join(lines)


def inject_ledger_block(raw_card_text: str, block: str, profile: dict) -> str:
    markers = profile.get("ledger_markers", {})
    begin_m = markers.get("begin", "<!-- LEDGER:BEGIN -->")
    end_m = markers.get("end", "<!-- LEDGER:END -->")

    pattern = re.compile(re.escape(begin_m) + r".*?" + re.escape(end_m), re.S)
    if pattern.search(raw_card_text):
        return pattern.sub(lambda _: block, raw_card_text, count=1)

    m1 = re.search(r"^## 1\b.*$", raw_card_text, re.M)
    if m1:
        m2 = re.search(r"^## (?!1\b)", raw_card_text[m1.end():], re.M)
        if m2:
            sec1_end = m1.end() + m2.start()
            sec1_body = raw_card_text[m1.end():sec1_end]
            sep_matches = list(re.finditer(r"^---\s*$", sec1_body, re.M))
            if sep_matches:
                insert_pos = m1.end() + sep_matches[-1].start()
                before = raw_card_text[:insert_pos].rstrip("\n")
                after = raw_card_text[insert_pos:]
                return before + "\n\n" + block + "\n\n" + after
            before = raw_card_text[:sec1_end].rstrip("\n")
            after = raw_card_text[sec1_end:]
            return before + "\n\n" + block + "\n\n" + after

    return raw_card_text.rstrip("\n") + "\n\n" + block + "\n"


def check_all(
    root: Path,
    card: Path,
    ledger_path: Path,
    profile: dict,
) -> list[str]:
    raw_card_text = _census.safe_read(card)
    expected_ledger, findings = scan_and_build_ledger(root, card, raw_card_text, profile)
    allowed_enum = set(profile.get("states", []))
    markers = profile.get("ledger_markers", {})
    begin_m = markers.get("begin", "<!-- LEDGER:BEGIN -->")
    end_m = markers.get("end", "<!-- LEDGER:END -->")

    if not ledger_path.is_file():
        findings.append(f"C10: LEDGER_FILE_MISSING: `{ledger_path}`")
        return findings

    try:
        actual_ledger = json.loads(ledger_path.read_text(encoding="utf-8-sig"))
    except Exception as exc:
        findings.append(f"C10: INVALID_LEDGER_JSON: {exc!r}")
        return findings

    actual_items = actual_ledger.get("items", [])
    if not isinstance(actual_items, list) or len(actual_items) == 0:
        findings.append("C9: EMPTY_LEDGER_ITEMS: `items` in status_ledger.json is empty")
        return findings

    for cf in CORE_FAMILIES:
        if sum(1 for it in actual_items if it.get("family") == cf) == 0:
            findings.append(f"C9: EMPTY_FAMILY_IN_LEDGER: `{cf}` has 0 items in status_ledger.json")

    for it in actual_items:
        p = str(it.get("path", ""))
        st = str(it.get("state", ""))
        fl = it.get("flags", []) or []
        ev = str(it.get("evidence", ""))
        prog = str(it.get("progress", ""))
        anchor = str(it.get("roadmap_anchor", ""))

        if st not in allowed_enum:
            findings.append(f"C6: OUT_OF_ENUM_STATE in ledger on `{p}`: state={st}")
        if (st == "TERMINATED") != ("archive" in Path(p).parts):
            findings.append(f"C7: TERMINATED_ARCHIVE_MISMATCH in ledger on `{p}`: state={st}")
        if st == "DELIVERED" and (not ev or "evidence_missing" in fl):
            findings.append(f"C8: DELIVERED_WITHOUT_EVIDENCE in ledger on `{p}`")
        if st == "PARTIAL":
            jm = re.match(r"^(\d+)\s*/\s*(\d+)", prog)
            if (jm and int(jm.group(1)) == 0) or (not prog and not ev):
                findings.append(f"C8: PARTIAL_WITHOUT_CLOSED_JALON in ledger on `{p}` (progress={prog})")
        if not anchor:
            findings.append(f"C8: MISSING_ROADMAP_ANCHOR in ledger on `{p}`")

    exp_map = {it["path"]: it for it in expected_ledger.get("items", [])}
    act_map = {str(it.get("path", "")): it for it in actual_items}
    if set(exp_map.keys()) != set(act_map.keys()):
        missing = sorted(set(exp_map.keys()) - set(act_map.keys()))
        extra = sorted(set(act_map.keys()) - set(exp_map.keys()))
        findings.append(
            f"C10: LEDGER_ITEM_SET_MISMATCH: missing={missing[:5]} extra={extra[:5]}"
        )

    for p, exp_it in sorted(exp_map.items()):
        act_it = act_map.get(p)
        if not act_it:
            continue
        for k in ("state", "progress", "evidence", "next_action", "flags", "roadmap_anchor", "family", "item_kind"):
            if act_it.get(k) != exp_it.get(k):
                findings.append(
                    f"C10: LEDGER_MISMATCH on `{p}` field `{k}`: ledger={act_it.get(k)!r} vs card/disk={exp_it.get(k)!r}"
                )

    expected_block = render_ledger_block(expected_ledger, profile)
    m_block = re.search(re.escape(begin_m) + r".*?" + re.escape(end_m), raw_card_text, re.S)
    if not m_block:
        findings.append("C10: MISSING_LEDGER_BLOCK in ROADMAP.md")
    elif m_block.group(0) != expected_block:
        findings.append("C10: ROADMAP_LEDGER_BLOCK_DRIFT in ROADMAP.md")

    return findings


def resolve_profile_path(root: Path, profile_arg: str) -> Path:
    p = Path(os.path.expanduser(profile_arg))
    if p.is_absolute() and p.is_file():
        return p
    cand1 = (root / p).resolve()
    if cand1.is_file():
        return cand1
    cand2 = (_REPO_ROOT / p).resolve()
    if cand2.is_file():
        return cand2
    cand3 = (
        _REPO_ROOT
        / "docs"
        / "agent"
        / "skills"
        / "scan-inventory"
        / "references"
        / "profile.example.json"
    ).resolve()
    return cand3


def main() -> int:
    ap = argparse.ArgumentParser(description="Status backbone ledger generator & C6-C10 guard.")
    ap.add_argument("--build", action="store_true", help="Generate status_ledger.json and update ROADMAP.md LEDGER block")
    ap.add_argument("--check", action="store_true", help="Validate C6-C10 guards against status_ledger.json and ROADMAP.md")
    ap.add_argument("--root", default=".")
    ap.add_argument("--card", default="docs/ROADMAP.md")
    ap.add_argument("--ledger", default="docs/agent/status_ledger.json")
    ap.add_argument(
        "--profile",
        default="docs/agent/skills/scan-inventory/references/profile.example.json",
    )
    a = ap.parse_args()

    root = Path(a.root).resolve()
    card = (root / a.card).resolve()
    ledger_path = (root / a.ledger).resolve()
    profile_path = resolve_profile_path(root, a.profile)

    if not profile_path.is_file():
        print("INSTRUMENT FAILURE: profile not found: " + str(profile_path), flush=True)
        return 2
    if not card.is_file():
        print("INSTRUMENT FAILURE: card not found: " + str(card), flush=True)
        return 2

    try:
        profile = json.loads(profile_path.read_text(encoding="utf-8-sig"))
    except Exception as exc:
        print("INSTRUMENT FAILURE: " + repr(exc), flush=True)
        return 2

    if a.build:
        try:
            raw_card_text = _census.safe_read(card)
            ledger, _ = scan_and_build_ledger(root, card, raw_card_text, profile)
            block = render_ledger_block(ledger, profile)
            new_card_text = inject_ledger_block(raw_card_text, block, profile)

            ledger_path.parent.mkdir(parents=True, exist_ok=True)
            json_bytes = (json.dumps(ledger, indent=2, sort_keys=True, ensure_ascii=False) + "\n").encode("utf-8")
            ledger_path.write_bytes(json_bytes)
            card.write_bytes(new_card_text.replace("\r\n", "\n").encode("utf-8"))
            print(
                f"BUILT: {ledger_path} ({ledger['summary']['total_items']} items) & updated {card}",
                flush=True,
            )
            if not a.check:
                return 0
        except Exception as exc:
            print("INSTRUMENT FAILURE: " + repr(exc), flush=True)
            return 2

    try:
        findings = check_all(root, card, ledger_path, profile)
    except Exception as exc:
        print("INSTRUMENT FAILURE: " + repr(exc), flush=True)
        return 2

    if findings:
        print(f"VERDICT: {len(findings)} C6-C10 finding(s):", flush=True)
        for f in findings:
            print("  - " + f, flush=True)
        return 1

    print("VERDICT: 0 finding(s) (C6-C10 GREEN)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

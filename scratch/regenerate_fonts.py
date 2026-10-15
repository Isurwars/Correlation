#!/usr/bin/env python3
"""
Regenerates RobotoData_ASCII_1.cpp, RobotoData_ASCII_2.cpp, and RobotoData_Extended.cpp
using static constexpr Point2D arrays and size_t stroke lengths.
"""

import re
import os

def split_args(arg_str):
    args = []
    current = []
    in_quote = False
    quote_char = ''
    brace_depth = 0
    i = 0
    while i < len(arg_str):
        ch = arg_str[i]
        if in_quote:
            current.append(ch)
            if ch == '\\':
                i += 1
                if i < len(arg_str):
                    current.append(arg_str[i])
            elif ch == quote_char:
                in_quote = False
        else:
            if ch in ("'", '"'):
                in_quote = True
                quote_char = ch
                current.append(ch)
            elif ch == '{':
                brace_depth += 1
                current.append(ch)
            elif ch == '}':
                brace_depth -= 1
                current.append(ch)
            elif ch == ',' and brace_depth == 0:
                args.append("".join(current).strip())
                current = []
            else:
                current.append(ch)
        i += 1
    if current:
        args.append("".join(current).strip())
    return args

def parse_glyphs(content):
    glyphs = []
    idx = 0
    while True:
        pos = content.find('add(', idx)
        if pos == -1:
            break
        paren_depth = 1
        i = pos + 4
        in_quote = False
        quote_char = ''
        while i < len(content) and paren_depth > 0:
            ch = content[i]
            if in_quote:
                if ch == '\\':
                    i += 1
                elif ch == quote_char:
                    in_quote = False
            else:
                if ch in ("'", '"'):
                    in_quote = True
                    quote_char = ch
                elif ch == '(':
                    paren_depth += 1
                elif ch == ')':
                    paren_depth -= 1
                    if paren_depth == 0:
                        break
            i += 1
        call_str = content[pos+4:i].strip()
        args = split_args(call_str)
        assert len(args) == 4, f"Expected 4 args in {call_str[:50]}, got {len(args)}: {args}"
        code, left, right, strokes_str = args
        
        strokes = []
        if strokes_str != '{}':
            inner = strokes_str
            if inner.startswith('{') and inner.endswith('}'):
                inner = inner[1:-1].strip()
            depth = 0
            start = -1
            for j, ch in enumerate(inner):
                if ch == '{':
                    if depth == 0:
                        start = j
                    depth += 1
                elif ch == '}':
                    depth -= 1
                    if depth == 0 and start != -1:
                        stroke_str = inner[start+1:j].strip()
                        pts = re.findall(r'\{\s*([^,]+)\s*,\s*([^}]+)\s*\}', stroke_str)
                        strokes.append(pts)
                        start = -1
        glyphs.append((code, left, right, strokes))
        idx = i + 1
    return glyphs

def code_to_id(code, index):
    if code.startswith("'") and code.endswith("'") and len(code) >= 3:
        inner = code[1:-1]
        if inner == '\\\\':
            return f"BACKSLASH_{index}"
        elif inner == "\\'":
            return f"SINGLEQUOTE_{index}"
        elif inner.isalnum():
            return f"{inner}_{index}"
        else:
            return f"SYM_{index}"
    elif code.startswith("0x"):
        return f"HEX_{code[2:]}_{index}"
    return f"GLYPH_{index}"

def generate_file(input_path, func_name, brief_desc):
    with open(input_path, 'r') as f:
        content = f.read()
    glyphs = parse_glyphs(content)
    
    filename = os.path.basename(input_path)
    lines = [
        "/**",
        f" * @file {filename}",
        f" * @brief {brief_desc}",
        " * @copyright Copyright © 2013-2026 Isaías Rodríguez (isurwars@gmail.com)",
        " * @par License",
        " * SPDX-License-Identifier: AGPL-3.0-only",
        " */",
        "",
        '#include "plotters/PathFont.hpp"',
        "",
        "namespace correlation::plotters {",
        "",
        "namespace {",
        ""
    ]
    
    glyph_info = []
    for idx, (code, left, right, strokes) in enumerate(glyphs):
        all_pts = []
        stroke_lens = []
        for s in strokes:
            stroke_lens.append(len(s))
            for pt in s:
                all_pts.append(pt)
        
        gid = code_to_id(code, idx)
        if all_pts:
            pts_id = f"GLYPH_PTS_{gid}"
            strokes_id = f"GLYPH_STROKES_{gid}"
            lines.append(f"constexpr Point2D {pts_id}[] = {{")
            chunk_size = 4
            for c_i in range(0, len(all_pts), chunk_size):
                chunk = all_pts[c_i:c_i+chunk_size]
                chunk_str = ", ".join(f"{{{pt[0]}, {pt[1]}}}" for pt in chunk)
                sep = "," if c_i + chunk_size < len(all_pts) else ""
                lines.append(f"    {chunk_str}{sep}")
            lines.append("};")
            lines.append("")
            stroke_lens_str = ", ".join(str(sl) for sl in stroke_lens)
            lines.append(f"constexpr size_t {strokes_id}[] = {{{stroke_lens_str}}};")
            lines.append("")
            glyph_info.append((code, left, right, pts_id, strokes_id))
        else:
            glyph_info.append((code, left, right, None, None))
            
    lines.append("} // namespace")
    lines.append("")
    lines.append(f"void Roboto::{func_name}() {{")
    for code, left, right, pts_id, strokes_id in glyph_info:
        if pts_id is None:
            lines.append(f"  add({code}, {left}, {right}, {{}}, {{}});")
        else:
            lines.append(f"  add({code}, {left}, {right}, {pts_id}, {strokes_id});")
    lines.append("}")
    lines.append("")
    lines.append("} // namespace correlation::plotters")
    lines.append("")
    
    with open(input_path, 'w') as f:
        f.write("\n".join(lines))
    print(f"Generated {input_path}: {len(glyphs)} glyphs.")

if __name__ == '__main__':
    generate_file('src/plotters/fonts/RobotoData_ASCII_1.cpp', 'registerBasic1', "Roboto font basic ASCII glyph data (Part 1: ' ' to 'L').")
    generate_file('src/plotters/fonts/RobotoData_ASCII_2.cpp', 'registerBasic2', "Roboto font basic ASCII glyph data (Part 2: 'M' to '~').")
    generate_file('src/plotters/fonts/RobotoData_Extended.cpp', 'registerExtended', "Roboto font extended Unicode glyph data.")

#!/usr/bin/env python3
import sys
import os
import traceback

def get_depth(name):
    stripped = name.lstrip()
    # Each indentation level is 2 spaces in standard Kraken reports
    return (len(name) - len(stripped)) // 2

def fix_report(report_file, contaminant_taxids, output_file):
    try:
        # Read the whole report
        lines = []
        if os.path.exists(report_file):
            with open(report_file, 'r', encoding='utf-8', errors='replace') as f:
                for line in f:
                    parts = line.rstrip('\r\n').split('\t')
                    if len(parts) < 6: continue
                    # 0: percentage, 1: clade_counts, 2: taxon_counts, 3: rank, 4: taxid, 5: name
                    lines.append(parts)

        if not lines:
            sys.stderr.write(f"Warning: No valid lines found in {report_file}. Creating empty output.\n")
            with open(output_file, 'w', encoding='utf-8') as f:
                pass
            return

        # 1. Zero out counts for identified contaminants
        for line in lines:
            if line[4] in contaminant_taxids:
                line[1] = '0' # clade counts
                line[2] = '0' # direct counts

        # 2. Build parent-child relationships
        child_indices = {}
        stack = [] # (depth, index)
        
        for i, line in enumerate(lines):
            name = line[5]
            depth = get_depth(name)
            
            while stack and stack[-1][0] >= depth:
                stack.pop()
            
            if stack:
                parent_idx = stack[-1][1]
                child_indices.setdefault(parent_idx, []).append(i)
            
            stack.append((depth, i))

        # 3. Iterative update of clade counts (post-order traversal-like)
        # We go backwards because the report is pre-order.
        # This way we process children before their parents.
        for i in range(len(lines) - 1, -1, -1):
            try:
                direct_count = int(lines[i][2])
            except ValueError:
                direct_count = 0
                
            clade_count = direct_count
            if i in child_indices:
                for child_idx in child_indices[i]:
                    try:
                        clade_count += int(lines[child_idx][1])
                    except ValueError:
                        pass
            lines[i][1] = str(clade_count)

        # 4. Total reads = sum of clade counts of all nodes at depth 0
        total_reads = 0
        for i, line in enumerate(lines):
            if get_depth(line[5]) == 0:
                try:
                    total_reads += int(line[1])
                except ValueError:
                    pass

        # 5. Update percentages
        if total_reads > 0:
            for line in lines:
                try:
                    line[0] = f"{int(line[1]) * 100 / total_reads:.5f}"
                except (ValueError, ZeroDivisionError):
                    line[0] = "0.00000"
        else:
            for line in lines:
                line[0] = "0.00000"

        # 6. Write the fixed report
        with open(output_file, 'w', encoding='utf-8') as f:
            for line in lines:
                f.write('\t'.join(line) + '\n')
                
    except Exception as e:
        sys.stderr.write(f"Error processing {report_file}: {str(e)}\n")
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: fix_kraken_report.py <report_file> <contaminants_file> <output_file>")
        sys.exit(1)
    
    report_file = sys.argv[1]
    contaminants_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Read contaminant taxids
    contam_taxids = set()
    if os.path.exists(contaminants_file):
        try:
            with open(contaminants_file, 'r', encoding='utf-8', errors='replace') as f:
                for line in f:
                    taxid = line.strip()
                    if taxid:
                        contam_taxids.add(taxid)
        except Exception as e:
            sys.stderr.write(f"Warning: Could not read contaminants file {contaminants_file}: {e}\n")
    
    fix_report(report_file, contam_taxids, output_file)

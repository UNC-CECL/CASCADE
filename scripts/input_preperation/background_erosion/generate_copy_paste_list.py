"""
Generate Copy-Paste Ready Background Erosion Rates List

This script reads your background_erosion_rates.csv and outputs a 
formatted Python list that you can copy directly into your hindcast script.

Output: A text file with the list ready to copy-paste
"""

import pandas as pd
import numpy as np

# ========================================================================
# CONFIGURATION
# ========================================================================

INPUT_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\data\background_erosion\background_erosion_rates.csv"
OUTPUT_FILE = r"C:\Users\hanna\PycharmProjects\CASCADE\data\background_erosion\1978_1997_copy_paste_list.txt"

# ========================================================================
# LOAD AND CONVERT
# ========================================================================

print("="*70)
print("BACKGROUND EROSION RATES → COPY-PASTE LIST GENERATOR")
print("="*70)

try:
    # Load the CSV
    df = pd.read_csv(INPUT_CSV)
    rates = df['background_rate_m_per_yr'].values
    
    print(f"\n✓ Loaded {len(rates)} background erosion rates")
    print(f"  Range: {rates.min():.3f} to {rates.max():.3f} m/yr")
    print(f"  Mean: {rates.mean():.4f} m/yr")
    
    # Convert to Python list (not NumPy array!)
    rates_list = rates.tolist()
    
    # Verify it's a Python list
    print(f"\n✓ Converted to Python list")
    print(f"  Type: {type(rates_list)}")
    print(f"  First element type: {type(rates_list[0])}")
    
    # ========================================================================
    # FORMAT FOR COPY-PASTE
    # ========================================================================
    
    # Create the formatted string
    output_lines = []
    output_lines.append("# Background erosion rates for CASCADE (m/yr)")
    output_lines.append("# Generated from: background_erosion_rates.csv")
    output_lines.append("# Domains 0-14: Left buffer (0.0)")
    output_lines.append("# Domains 15-104: Real island (calculated rates)")
    output_lines.append("# Domains 105-119: Right buffer (0.0)")
    output_lines.append("")
    output_lines.append("BACKGROUND_EROSION_RATES = [")
    
    # Format rates with 10 values per line for readability
    for i in range(0, len(rates_list), 10):
        # Get next 10 values (or remaining values)
        chunk = rates_list[i:i+10]
        
        # Format each value to 4 decimal places
        formatted_chunk = [f"{val:7.4f}" for val in chunk]
        
        # Create line with comments showing domain indices
        line = "    " + ", ".join(formatted_chunk) + ","
        
        # Add comment showing which domains
        if i == 0:
            line += "  # Domains 0-9"
        elif i + 10 >= len(rates_list):
            line += f"  # Domains {i}-{len(rates_list)-1}"
        else:
            line += f"  # Domains {i}-{i+9}"
        
        output_lines.append(line)
    
    output_lines.append("]")
    output_lines.append("")
    output_lines.append("# To use: Replace the line in your hindcast script:")
    output_lines.append("# BACKGROUND_EROSION_RATES = [0.0] * TOTAL_DOMAINS")
    output_lines.append("# with the list above (including the BACKGROUND_EROSION_RATES = [ ... ] part)")
    
    # ========================================================================
    # SAVE TO FILE
    # ========================================================================
    
    output_text = "\n".join(output_lines)
    
    with open(OUTPUT_FILE, 'w') as f:
        f.write(output_text)
    
    print(f"\n✓ Saved formatted list to: {OUTPUT_FILE}")
    
    # ========================================================================
    # ALSO PRINT TO CONSOLE (for immediate copy-paste)
    # ========================================================================
    
    print("\n" + "="*70)
    print("COPY-PASTE READY LIST (also saved to file):")
    print("="*70)
    print()
    print(output_text)
    print()
    print("="*70)
    
    # ========================================================================
    # STATISTICS
    # ========================================================================
    
    print("\nVERIFICATION:")
    print(f"  Total values: {len(rates_list)}")
    print(f"  Buffer domains (0-14): {sum(1 for r in rates_list[0:15] if r == 0.0)} zeros")
    print(f"  Real island (15-104): Range {min(rates_list[15:105]):.3f} to {max(rates_list[15:105]):.3f} m/yr")
    print(f"  Buffer domains (105-119): {sum(1 for r in rates_list[105:120] if r == 0.0)} zeros")
    
    # Check if any values are problematic
    if any(abs(r) > 10 for r in rates_list):
        print("\n⚠️  WARNING: Some rates exceed ±10 m/yr - check if this is expected")
    else:
        print("\n✓ All rates within reasonable range (±10 m/yr)")
    
    print("\n" + "="*70)
    print("DONE!")
    print("="*70)
    print("\nNext steps:")
    print(f"1. Open: {OUTPUT_FILE}")
    print("2. Copy everything between (and including):")
    print("   BACKGROUND_EROSION_RATES = [")
    print("   ...")
    print("   ]")
    print("3. Paste into your hindcast script, replacing:")
    print("   BACKGROUND_EROSION_RATES = [0.0] * TOTAL_DOMAINS")
    print("="*70)

except FileNotFoundError:
    print(f"\n❌ ERROR: File not found: {INPUT_CSV}")
    print("   Check the path is correct.")
except KeyError:
    print(f"\n❌ ERROR: Column 'background_rate_m_per_yr' not found in CSV")
    print("   Check the CSV file has the correct column.")
except Exception as e:
    print(f"\n❌ ERROR: {e}")
    import traceback
    traceback.print_exc()

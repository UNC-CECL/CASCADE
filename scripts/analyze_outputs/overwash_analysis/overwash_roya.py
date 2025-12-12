import csv

output_rows = []
output_file = "overwash_flux_by_domain_year.csv"

# Loop through years and domains
for year in range(0, 20):
    for domain_id, domain_segment in enumerate(cascade.barrier3d):

        # Limit to real Pea Island domain range (1 to 71)
        if 1 <= domain_id <= 71:

            # Ensure the record exists for this year
            if year < len(domain_segment._QowTS):
                overwash_flux = domain_segment._QowTS[year]
            else:
                overwash_flux = None

            # Store results
            output_rows.append({
                "domain": domain_id,
                "year": year,
                "overwash_flux": overwash_flux
            })

# Write CSV
with open(output_file, "w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=["domain", "year", "overwash_flux"])
    writer.writeheader()
    writer.writerows(output_rows)

print(f"\nCSV successfully written → {output_file}")
# Nourishment datasets

- `NC_beachno-episodes-2025-10-28.xlsx` — the North Carolina episodes from the national beach nourishment database, downloaded 2025-10-28. Left as downloaded.
- `Hatteras_BN_data.xlsx` — the Hatteras Island rows taken from it, plus corrections.

## Corrections to the source dataset

**2017–18 Buxton nourishment is missing.** The source has no record of it. It placed 2.6 million cubic yards along 2.94 miles of northern Buxton, from the Haulover Day Use Area south to the groin at the old Cape Hatteras Lighthouse site. Work began in summer 2017 (about 46% placed by November 2017) and finished on 27 February 2018 (Weeks Marine, about $22 million). It was added to `Hatteras_BN_data.xlsx` on 2026-10-02 with `yearCompleted` 2018. The hindcast places it in 2017 (`HATTERAS_NOURISHMENT_PROJECTS`).

**The source's 2022 Buxton row (id 4274) cites the wrong project.** Its volume (1.2 million cubic yards) and year match the 2022 renourishment, but its source URL is the March 2018 article reporting the end of the 2017–18 project. That row is now in `Hatteras_BN_data.xlsx` with its values unchanged and its source replaced. The original URL is recorded in `otherInfo`.

## Additions not yet in the national database

- **Avon 2026:** about 375,000 cubic yards along about 1 mile, from just south of Avon Pier to the south village limit, placed 27 May to 25 June 2026.
- **Buxton 2026:** 2.0 million cubic yards **planned**, Haulover to the lighthouse groin field. Pumping began 31 July 2026, with about 75% placed by 16 September. Replace with the as-built volume once it's published.

Neither project's cost has been reported separately (the combined contract is about $45 million), so their cost fields are 0, the sheet's convention for an unknown cost.

Fills ruled out on 2026-10-02 (not separate placements): FEMA reimbursements after Florence and Dorian for Buxton paid for sand inside the 2022 project; state highway (NCDOT) dune pushing and sandbags at the S-curves, Rodanthe and north Buxton moved island sand, not imported fill.

Sources:
- https://coastalreview.org/2018/03/buxton-beach-nourishment-project-complete/
- https://outerbanksvoice.com/2018/03/01/delayed-buxton-beach-nourishment-project-is-finally-done/
- https://www.outerbanksvoice.com/2017/11/06/buxton-beach-project-will-miss-contract-deadline-by-2-months/
- https://content.govdelivery.com/accounts/NCDARECOUNTY/bulletins/3285a2f (2022 renourishment)
- https://www.outerbanksvoice.com/2026/05/26/avon-beach-nourishment-project-expected-to-begin-on-wednesday-may-27/ (Avon 2026)
- https://www.outerbanksvoice.com/2026/09/17/latest-update-on-buxton-beach-nourishment-and-groin-repair/ (Buxton 2026)
- https://islandfreepress.org/blog/avon-and-buxton-beach-nourishment-faqs-2026-edition/ (2026 projects)

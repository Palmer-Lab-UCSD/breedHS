# breedHS code instances

Use these instances as templates for running breedHS. They are organized by institution and pairing protocol:
- `hsw`: Use this code for in-house pairing at HS West (using only rats produced at HSW)
- `wfu`: Use this code for in-house pairing at Wake Forest (using only rats produced at WFU)
- `hsw_exchange`: Use this code for pairing exchange generations at HS West (when HSW receives breeders from WFU)
- `wfu_exchange`: Use this code for pairing exchange generations at Wake Forest (when WFU receives breeders from HSW)

## HS West
Code in `hsw` uploaded 1/21/26. This code should work as-is (after filling out the arguments file).

## Wake Forest
 No code present in `wfu` yet. Some functions have been written and I began testing hypothetical WFU pairings during HSW gen104, but status is incomplete. This code base should be finished ASAP

## HSW Exchange
Code in `hsw_exchange` should be stable and usable as-is after filling out the arguments file. Uploaded 1/21/26 following HSW gen108_109 pairings.

This code base pairs HSW rats during exchanges, from the perspective of HS West. Use this code only when HS West receives breeders from Wake Forest. Note that script `03_propose_pairs` lacks steps for identifying alternative/replacement pairs, due to complications imposed by WFU IDs. 

## WFU Exchange
Code uploaded 1/21/26 using code from WFU gen51_52 pairings. This might not be 100% usable as-is.

WFU exchanges require an early special step, `02_fams_to_ship`, which identifies best-candidate families to send from HSW to WFU based on shared kinship with the **previous** WFU generation. This is because HSW animals must be shipped to WFU before all WFU rats are born, thus current kinship can't be estimated. The previous generation's kinship is used as a proxy for prioritizing low-kinship HSW rats to ship.

During WFU exchanges, execute `01_make_id_maps` and `02_fams_to_ship` to identify breeders to provide from HSW. When beginning official pairing, re-run `01_make_id_maps`, then proceed to `03_setup_peds`. 

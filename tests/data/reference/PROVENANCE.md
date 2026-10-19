# Reference outputs for the integration tests (provenance)

These files are byte-identical copies of the scientific outputs of the
Phase 2D audit baseline of LorenzCycleToolkit 2.0.0 (commit
d38cda7e37d8e8a3a937a5919640a94bef19e34a). They are used by
`tests/output_checks.py` to compare the integration-test outputs value by
value. They are a reproducibility reference for this input and this
configuration, not an independent physical validation.

| Directory | Audit run | Command (from the repository root) |
|---|---|---|
| `ERA5_fixed/` | `run_2d/fixed_run1` | `lorenzcycletoolkit.py samples/testdata_ERA5.nc -f -r --box_limits inputs/box_limits_Reg1` |
| `ERA5_track/` | `run_2d/track_original_run1` | `lorenzcycletoolkit.py samples/testdata_ERA5.nc -t -r --trackfile inputs/track_testdata_ERA5` |

Inputs: `samples/testdata_ERA5.nc` (sha256 54347c9d93041ab4cc686305580893f2f510f8896f7cf23baf8718f933542005),
`inputs/namelist_ERA5` (575825542c7adf6a3a98f3a43d39811fa26d89a3daaed9e4e6c9faa2aea9bc4f),
`inputs/box_limits_Reg1` (dd301b7a98b35d055823e3121cbb37d246343e77be1a436c4e62ee5dbe135eea),
`inputs/track_testdata_ERA5` (5b82ff3cd827b41ee670e3625bd3486fdf6b8edfc4fad1e095844eb50e435358).

Environment of the audit runs: Python 3.12.15, x86_64 (Rosetta) on macOS 14.1.2, numpy 2.0.0,
scipy 1.14.0, pandas 2.2.2, xarray 2024.2.0, netCDF4 1.7.1 (netCDF-C 4.9.2, HDF5 1.14.3),
metpy 1.6.2, pint 0.24.3. Both audit repetitions were bit-identical.

Tolerances applied by the comparison are those registered for the audit
(|a - b| <= atol + rtol * max(|a|, |b|), rtol = 1e-10): energies 1e-6 J m-2;
rates 1e-9 W m-2; column M 1e-10 kg m-2 s-1; per-level M and vorticity
1e-12 s-1; coordinates 1e-9 degree with rtol = 0; other values 1e-12.
Times, level headers and text cells must be identical.

`MANIFEST.tsv` lists every reference file with its case, class, units,
relation to the column integral, source run and path in the audit, sha256,
code commit and command. These 49 files are regression fixtures for the
tests; the complete audit baseline (154 files, including the moving-track
case and the run logs) stays in the audit area and is not replaced.

No approved numerical reference exists for the NCEP-R2 cases; their
outputs are checked structurally only.

sha256 of every file:

```
6f203e26ee8550917b1a0d2e42387fa685dd18fdfbb070410b5b6a908255e8a7  ./ERA5_fixed/results_vertical_levels/Ae_level.csv
f9047f74dca54a3a8e844a058ec6b434b2e1598791f00b5ee3a81186d3cef728  ./ERA5_fixed/results_vertical_levels/Az_level.csv
9a5e48def3bb9120df476e08a1ae46b36fc4820f9ed6c6bbba7d3dc684489164  ./ERA5_fixed/results_vertical_levels/C_overturning_level.csv
cbcd33ea836d9fdad48dcfe4041289ac66fc151643db7cf85905ffbf88efe0a7  ./ERA5_fixed/results_vertical_levels/Ca_1_level.csv
4552e635c53ba306ef502c937befa202b2ab9bfe70ac93c5aca1fedc5b969c5e  ./ERA5_fixed/results_vertical_levels/Ca_2_level.csv
2e06476b850f54a84ef1c8300909f816cd3693bcf9bf254f4546be8fdd40db3f  ./ERA5_fixed/results_vertical_levels/Ca_level.csv
46a775644107ce5df5be9f4ad763e4d74a1dd5ab757a0c9631ff7c8237e0624c  ./ERA5_fixed/results_vertical_levels/Ce_1_level.csv
05440abb184c7d9dae0684a7d41fe407be5baa9213cae25df191c440b2731b37  ./ERA5_fixed/results_vertical_levels/Ce_2_level.csv
204f41b148f6a6370be1af769aa3feac4a564f49f2c6cb3c19a29c10a32151de  ./ERA5_fixed/results_vertical_levels/Ce_level.csv
7e382849a432ec459570917a1526b656f9869169657f3144410381ec3bc652bb  ./ERA5_fixed/results_vertical_levels/Ck_1_level.csv
3871d023a86ce892b3058fa0586dd63dd4bea85ef7610973c0dd60ace4377782  ./ERA5_fixed/results_vertical_levels/Ck_2_level.csv
eb5e346c607e99b048efbfeeca8162f482a6b4e9192720ff6b420fff8cd45b98  ./ERA5_fixed/results_vertical_levels/Ck_3_level.csv
d4d5968c1c8fdc0835ab6484efcdf27cf62a6a9bd1902a06b4a4f6d41e81acca  ./ERA5_fixed/results_vertical_levels/Ck_4_level.csv
f0598a7888334d9520c9519c5c4f34197e93e2ac715d40db572e28e99ab4fdb6  ./ERA5_fixed/results_vertical_levels/Ck_5_level.csv
2b37ea2f00083c1998cb068c4dc839781e7e8365de332ea2e8a569b46abb5d3b  ./ERA5_fixed/results_vertical_levels/Ck_level.csv
aff4304aff285e0c850282a43063e815aa8fbe3138cea762df059acfd094edb7  ./ERA5_fixed/results_vertical_levels/Cz_1_level.csv
853343abfc3e4a1707e154fa37e31dc9444a68204917f53d2982007cea79f117  ./ERA5_fixed/results_vertical_levels/Cz_2_level.csv
a3b4b4e6190a3a7c2f853ed7bda2f336ad62896e4b5f10602dc00a6b244652ec  ./ERA5_fixed/results_vertical_levels/Cz_level.csv
6126d3da8050d1e02561b62e88628f5dc25c9763a2e96f3b248e589b9d77009e  ./ERA5_fixed/results_vertical_levels/Ge_level.csv
4770deb27d80f88284aac780dd5a6481d10c44f9f344fe5a0c7cc89488775a16  ./ERA5_fixed/results_vertical_levels/Gz_level.csv
94da27b3ce47b25d0121d84798ccb9b1b15795208754229744266994b8726e1a  ./ERA5_fixed/results_vertical_levels/Ke_level.csv
2d2ff355c8caab9bfc2a7d11c94bf5f64412be71cb95feb1f19fbd81b5148445  ./ERA5_fixed/results_vertical_levels/Kz_level.csv
422ae7f27188639775bf89fb926860f977e3773d4e1ea023ffb72142242e0e1b  ./ERA5_fixed/results_vertical_levels/M_level.csv
5e5f1068a582b7e898c39d2e5dcee21c91879271f744d5b54de5969e20da0989  ./ERA5_fixed/testdata_ERA5_fixed_results.csv
f0cc1fc5cbb1725b20c994d087a531645712e9ce0a71b0bbd5791031b9c7cf94  ./ERA5_track/results_vertical_levels/Ae_level.csv
3ec9ae4c4dc485507852ce4c7ddbd105aba4f5e5d994eedad272b527e47e92c9  ./ERA5_track/results_vertical_levels/Az_level.csv
9cc8db0d103a0004c65a209172a686a56ce798e597bd99a314a3b1d319810341  ./ERA5_track/results_vertical_levels/C_overturning_level.csv
29c6b441b3f92abda1a455270401fb2e40eb49dede8b5b51e2daecd573e60d01  ./ERA5_track/results_vertical_levels/Ca_1_level.csv
cbd6c217821a2f9bc3948d7ee635d3f51b538341c84927d88fcaad22f91cdd38  ./ERA5_track/results_vertical_levels/Ca_2_level.csv
b44b54984922e2f7194057b869ab7b7d197563bb78233ebc60ba5b3e70ddffee  ./ERA5_track/results_vertical_levels/Ca_level.csv
60d7d744a55f75e23a2621271fd7e7b3c1fd9ba9a4f3a50a027b580b1fa8e28f  ./ERA5_track/results_vertical_levels/Ce_1_level.csv
4439f9f412e71f3b611cea2ab12d52ca232d962af8a2ac36c4dad055d2e000c0  ./ERA5_track/results_vertical_levels/Ce_2_level.csv
bbdc5e7a53e2eb8997ffc114fe3899f0e041dd7eb39525f25e1d9c7f40452d4c  ./ERA5_track/results_vertical_levels/Ce_level.csv
248e598073e0ac5677e2eac4dd915ff6e58cbf29a3831a26078cfec90bd45111  ./ERA5_track/results_vertical_levels/Ck_1_level.csv
450af2f97386823a5289f8be70e14021eb177b86a0092a956471132e2f424938  ./ERA5_track/results_vertical_levels/Ck_2_level.csv
8aee690220735308f71a0265d2e102fe2037edfcb93196e9134a67563fa4efff  ./ERA5_track/results_vertical_levels/Ck_3_level.csv
89d7b6431d85b99c41a6d40bb595ec0d25afab159c3563971520877a3c5cec0b  ./ERA5_track/results_vertical_levels/Ck_4_level.csv
87ccc8b89fb19cd96271117b12f80abc537954808fcb23ff7308ec3cd8f8e843  ./ERA5_track/results_vertical_levels/Ck_5_level.csv
3ce9640c434bc40a8109e230ab81a20098ed354805491bd9a25747c18a26c9d1  ./ERA5_track/results_vertical_levels/Ck_level.csv
60d7d744a55f75e23a2621271fd7e7b3c1fd9ba9a4f3a50a027b580b1fa8e28f  ./ERA5_track/results_vertical_levels/Cz_1_level.csv
56656a26c44d5f6e4ff3569fcd2585e6c809ed959d0c5dd45b9c9707ad16bfe6  ./ERA5_track/results_vertical_levels/Cz_2_level.csv
61907878d495f1d6096deb3b0a0bb36b4124f2a5cb190de24df4c09ba97eb15f  ./ERA5_track/results_vertical_levels/Cz_level.csv
7e2a642da1257d277bc64acd5e2aa82e61b91ba5907f12a67edf2556624d7cdf  ./ERA5_track/results_vertical_levels/Ge_level.csv
1e8faede9c39f30839931c5bf6b2a8bd0e119708d3e85fa5ec06e2e0ef347453  ./ERA5_track/results_vertical_levels/Gz_level.csv
0266fb45bc228b4af011ef3673fd81d54c4f4d58f50043d6b0c4ffee92887d1b  ./ERA5_track/results_vertical_levels/Ke_level.csv
0c3384b5a2d0b23660ba2d49ca03ed39898950cf4e5baa82c3ad8ea0ea869e62  ./ERA5_track/results_vertical_levels/Kz_level.csv
4cb7614386cc8dfdf1ee646184368d7450306041b6348c2136756ad906a9fe3b  ./ERA5_track/results_vertical_levels/M_level.csv
a5144092f5122f2211716df18e9f8d642edd7aa3c3f769c0914fd6ee8b17be49  ./ERA5_track/testdata_ERA5_track_results.csv
e8cdc6066b89f9bc251d58b2a7abcd101c1d1194cc7f35216eceed5b0dd8c0ec  ./ERA5_track/testdata_ERA5_track_trackfile
```

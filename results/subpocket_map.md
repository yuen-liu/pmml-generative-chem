# KLIFS position -> KinFragLib subpocket lookup table

Subpocket centres: mean C-alpha of KinFragLib anchor residues (AP: [46, 51, 75, 15], FP: [72, 51, 10, 81], SE: [51], GA: [45, 17, 81], B1: [81, 28, 43, 38], B2: [18, 24, 70, 83]), from Sydow et al., *J. Chem. Inf. Model.* 2020, 60, 6081-6094 (doi:10.1021/acs.jcim.0c00839), code at github.com/volkamerlab/KinaseFocusedFragmentLibrary @88b61d8 (`fragmentation/cli.py` L28-33). Each residue is assigned to the centre nearest to any of its heavy atoms, or to `other` if its closest heavy atom is >= 8.0 A from all centres (KinFragLib's X-pool cutoff). The table gives the modal label over 3619 main-set structures; `modal_share` is that label's frequency. Segment names are KLIFS's (van Linden et al., *J. Med. Chem.* 2014) and are shown for reference only.

|   klifs_pos | klifs_segment   | subpocket   |   modal_share | runner_up   |   runner_up_share |   n_structures |   median_d_nearest_A | note                          |
|------------:|:----------------|:------------|--------------:|:------------|------------------:|---------------:|---------------------:|:------------------------------|
|           1 | I               | other       |         0.96  | SE          |             0.04  |           3614 |                10.99 |                               |
|           2 | I               | other       |         0.999 | SE          |             0.001 |           3617 |                10.92 |                               |
|           3 | I               | AP          |         0.961 | other       |             0.017 |           3618 |                 6.1  |                               |
|           4 | g.l             | other       |         0.521 | FP          |             0.47  |           3579 |                 8.05 | unstable (modal share < 0.75) |
|           5 | g.l             | FP          |         0.509 | other       |             0.483 |           3480 |                 7.93 | unstable (modal share < 0.75) |
|           6 | g.l             | other       |         0.594 | FP          |             0.374 |           3351 |                 8.44 | unstable (modal share < 0.75) |
|           7 | g.l             | other       |         0.571 | B2          |             0.295 |           3247 |                 8.56 | unstable (modal share < 0.75) |
|           8 | g.l             | B2          |         0.636 | other       |             0.229 |           3291 |                 6.12 | unstable (modal share < 0.75) |
|           9 | g.l             | other       |         0.664 | B2          |             0.269 |           3442 |                 8.66 | unstable (modal share < 0.75) |
|          10 | II              | GA          |         0.701 | other       |             0.297 |           3592 |                 7.69 | unstable (modal share < 0.75) |
|          11 | II              | GA          |         0.865 | AP          |             0.094 |           3618 |                 5.58 |                               |
|          12 | II              | other       |         0.591 | GA          |             0.296 |           3619 |                 8.08 | unstable (modal share < 0.75) |
|          13 | II              | other       |         0.976 | AP          |             0.014 |           3619 |                 8.99 |                               |
|          14 | III             | AP          |         0.802 | other       |             0.192 |           3619 |                 7.74 |                               |
|          15 | III             | AP          |         0.83  | GA          |             0.17  |           3619 |                 4.88 |                               |
|          16 | III             | GA          |         1     |             |             0     |           3619 |                 5.77 |                               |
|          17 | III             | GA          |         0.603 | B2          |             0.397 |           3618 |                 3.31 | unstable (modal share < 0.75) |
|          18 | III             | GA          |         0.969 | other       |             0.031 |           3614 |                 7.24 |                               |
|          19 | III             | B2          |         0.655 | other       |             0.248 |           3588 |                 7.01 | unstable (modal share < 0.75) |
|          20 | αC              | B2          |         0.728 | other       |             0.272 |           3576 |                 6.19 | unstable (modal share < 0.75) |
|          21 | αC              | B2          |         0.8   | other       |             0.183 |           3597 |                 6.44 |                               |
|          22 | αC              | other       |         0.682 | B2          |             0.295 |           3603 |                 8.44 | unstable (modal share < 0.75) |
|          23 | αC              | B2          |         0.752 | other       |             0.247 |           3609 |                 7.26 |                               |
|          24 | αC              | B2          |         0.938 | other       |             0.031 |           3612 |                 2.41 |                               |
|          25 | αC              | B1          |         0.906 | B2          |             0.092 |           3619 |                 5.68 |                               |
|          26 | αC              | B1          |         0.767 | other       |             0.233 |           3616 |                 7.68 |                               |
|          27 | αC              | B1          |         0.831 | other       |             0.137 |           3619 |                 6.95 |                               |
|          28 | αC              | B1          |         0.951 | B2          |             0.048 |           3619 |                 1.84 |                               |
|          29 | αC              | B1          |         0.995 | other       |             0.005 |           3619 |                 6.24 |                               |
|          30 | αC              | other       |         0.931 | B1          |             0.069 |           3615 |                 8.74 |                               |
|          31 | b.l             | B1          |         0.528 | other       |             0.472 |           3616 |                 7.96 | unstable (modal share < 0.75) |
|          32 | b.l             | other       |         0.813 | B1          |             0.187 |           3617 |                10.06 |                               |
|          33 | b.l             | other       |         0.999 | FP          |             0.001 |           3598 |                11.56 |                               |
|          34 | b.l             | other       |         0.999 | FP          |             0.001 |           3619 |                10.81 |                               |
|          35 | b.l             | B1          |         0.639 | other       |             0.36  |           3617 |                 7.81 | unstable (modal share < 0.75) |
|          36 | b.l             | AP          |         0.71  | B1          |             0.237 |           3619 |                 4.77 | unstable (modal share < 0.75) |
|          37 | b.l             | B1          |         0.998 | FP          |             0.001 |           3619 |                 5.41 |                               |
|          38 | IV              | B1          |         1     |             |             0     |           3619 |                 3.48 |                               |
|          39 | IV              | B1          |         0.999 | FP          |             0.001 |           3619 |                 6.36 |                               |
|          40 | IV              | B1          |         0.999 | other       |             0.001 |           3619 |                 6.8  |                               |
|          41 | IV              | B1          |         0.824 | other       |             0.176 |           3619 |                 7.45 |                               |
|          42 | V               | GA          |         0.711 | other       |             0.287 |           3618 |                 7.8  | unstable (modal share < 0.75) |
|          43 | V               | B1          |         0.765 | GA          |             0.235 |           3618 |                 4.02 |                               |
|          44 | V               | B1          |         0.58  | GA          |             0.42  |           3619 |                 5.7  | unstable (modal share < 0.75) |
|          45 | GK              | GA          |         0.922 | B1          |             0.076 |           3619 |                 1.49 |                               |
|          46 | hinge           | AP          |         0.999 | FP          |             0.001 |           3618 |                 4.3  |                               |
|          47 | hinge           | AP          |         0.999 | FP          |             0.001 |           3619 |                 4.91 |                               |
|          48 | hinge           | AP          |         0.874 | SE          |             0.125 |           3619 |                 3.59 |                               |
|          49 | linker          | SE          |         0.952 | other       |             0.04  |           3606 |                 4.46 |                               |
|          50 | linker          | SE          |         1     |             |             0     |           2295 |                 2.43 |                               |
|          51 | linker          | SE          |         1     |             |             0     |           3608 |                 0    |                               |
|          52 | linker          | SE          |         1     | FP          |             0     |           3613 |                 2.42 |                               |
|          53 | αD              | SE          |         0.999 | FP          |             0.001 |           3619 |                 5.65 |                               |
|          54 | αD              | SE          |         0.939 | FP          |             0.054 |           3616 |                 7.57 |                               |
|          55 | αD              | SE          |         0.948 | other       |             0.051 |           3618 |                 5.73 |                               |
|          56 | αD              | SE          |         0.952 | other       |             0.048 |           3617 |                 6.26 |                               |
|          57 | αD              | other       |         0.99  | SE          |             0.01  |           3476 |                 9.13 |                               |
|          58 | αD              | other       |         0.999 | SE          |             0.001 |           3471 |                10.73 |                               |
|          59 | αD              | other       |         0.81  | SE          |             0.19  |           3600 |                 9.48 |                               |
|          60 | αE              | other       |         1     |             |             0     |           3614 |                12.3  |                               |
|          61 | αE              | other       |         0.999 | B2          |             0.001 |           3614 |                11.02 |                               |
|          62 | αE              | other       |         1     |             |             0     |           3615 |                15.74 |                               |
|          63 | αE              | other       |         1     |             |             0     |           3614 |                17.45 |                               |
|          64 | αE              | other       |         1     |             |             0     |           3610 |                14    |                               |
|          65 | VI              | other       |         1     |             |             0     |           3614 |                15.14 |                               |
|          66 | VI              | other       |         0.933 | B2          |             0.067 |           3614 |                10.76 |                               |
|          67 | VI              | other       |         0.93  | B2          |             0.07  |           3614 |                10.19 |                               |
|          68 | c.l             | B2          |         0.701 | other       |             0.25  |           3615 |                 7.77 | unstable (modal share < 0.75) |
|          69 | c.l             | other       |         0.968 | B2          |             0.031 |           3616 |                 9.98 |                               |
|          70 | c.l             | other       |         0.429 | B2          |             0.42  |           3616 |                 7.91 | unstable (modal share < 0.75) |
|          71 | c.l             | other       |         0.999 | FP          |             0.001 |           3619 |                10.02 |                               |
|          72 | c.l             | FP          |         0.902 | other       |             0.098 |           3619 |                 7.2  |                               |
|          73 | c.l             | FP          |         0.989 | SE          |             0.007 |           3619 |                 6.78 |                               |
|          74 | c.l             | FP          |         0.999 | SE          |             0.001 |           3618 |                 3.1  |                               |
|          75 | c.l             | FP          |         0.999 | AP          |             0.001 |           3619 |                 3.61 |                               |
|          76 | VII             | FP          |         0.96  | SE          |             0.04  |           3619 |                 5.64 |                               |
|          77 | VII             | AP          |         0.947 | FP          |             0.046 |           3617 |                 2.46 |                               |
|          78 | VII             | SE          |         0.999 | AP          |             0.001 |           3619 |                 4.66 |                               |
|          79 | VIII            | FP          |         0.518 | AP          |             0.24  |           3619 |                 7.69 | unstable (modal share < 0.75) |
|          80 | xDFG            | FP          |         0.933 | GA          |             0.042 |           3619 |                 4.08 |                               |
|          81 | xDFG            | B2          |         0.718 | FP          |             0.239 |           3619 |                 3.2  | unstable (modal share < 0.75) |
|          82 | xDFG            | B2          |         0.993 | B1          |             0.007 |           3619 |                 3.44 |                               |
|          83 | xDFG            | B2          |         0.999 | other       |             0.001 |           3551 |                 2.46 |                               |
|          84 | a.l             | B2          |         0.997 | FP          |             0.001 |           3451 |                 4.33 |                               |
|          85 | a.l             | B2          |         0.969 | other       |             0.025 |           3337 |                 6.24 |                               |

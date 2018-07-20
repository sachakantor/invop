MAX
    700 prem_pet_hn + 700 prem_pet_mn + 700 prem_pet_ln + 700 prem_pet_rg + 700 prem_pet_cg + 600 reg_pet_hn + 600 reg_pet_mn + 600 reg_pet_ln + 600 reg_pet_rg + 600 reg_pet_cg + 400 jf_ho + 400 jf_lo + 400 jf_co + 400 jf_res + 350 fo_ho + 350 fo_lo + 350 fo_co + 350 fo_res + 150 lube_res
ST.
    d1 + d2 <= 45000
    rg_ln + rg_mn + rg_hn <= 10000
    cracked_lo + cracked_ho <= 8000
    0.4 reg_pet_ln + 0.4 reg_pet_mn + 0.4 reg_pet_hn + 0.4 reg_pet_cg + 0.4 reg_pet_rg - prem_pet_ln - prem_pet_mn - prem_pet_hn - prem_pet_rg - prem_pet_cg <= 0
    reg_pet_ln + prem_pet_ln + rg_ln - 0.1 d1 - 0.15 d2 <= 0
    reg_pet_mn + prem_pet_mn + rg_mn - 0.2 d1 - 0.25 d2 <= 0
    reg_pet_hn + prem_pet_hn + rg_hn - 0.2 d1 - 0.18 d2 <= 0
    cracked_lo + jf_lo + fo_lo - 0.12 d1 - 0.08 d2 <= 0
    cracked_ho + jf_ho + fo_ho - 0.2 d1 - 0.19 d2 <= 0
    jf_res + fo_res + lube_res - 0.13 d1 - 0.12 d2 <= 0
    reg_pet_rg + prem_pet_rg - 0.6 rg_ln - 0.52 rg_mn - 0.45 rg_hn <= 0
    reg_pet_cg + prem_pet_cg - 0.28 cracked_lo - 0.2 cracked_ho <= 0
    fo_co + jf_co - 0.68 cracked_lo - 0.75 cracked_ho <= 0
    jf_lo + 0.6 jf_ho + 1.5 jf_co + 0.05 jf_res - jf_lo - jf_ho - jf_co - jf_res <= 0
    fo_lo - 10 fo_res = 0
    fo_co - 4 fo_res = 0
    fo_ho - 3 fo_res = 0
    70 reg_pet_hn + 80 reg_pet_mn + 90 reg_pet_ln + 115 reg_pet_rg + 105 reg_pet_cg - 84 reg_pet_hn - 84 reg_pet_mn - 84 reg_pet_ln - 84 reg_pet_rg - 84 reg_pet_cg => 0
    70 prem_pet_hn + 80 prem_pet_mn + 90 prem_pet_ln + 115 prem_pet_rg + 105 prem_pet_cg - 94 prem_pet_hn - 94 prem_pet_mn - 94 prem_pet_ln - 94 prem_pet_rg - 94 prem_pet_cg => 0
BOUNDS
    0 <= d1 <= 20000
    0 <= d2 <= 30000
    1000 <= lube_res <= 2000
    0 <= prem_pet_hn
    0 <= prem_pet_mn
    0 <= prem_pet_ln
    0 <= prem_pet_rg
    0 <= prem_pet_cg
    0 <= reg_pet_hn
    0 <= reg_pet_mn
    0 <= reg_pet_ln
    0 <= reg_pet_rg
    0 <= reg_pet_cg
    0 <= jf_ho
    0 <= jf_lo
    0 <= jf_co
    0 <= jf_res
    0 <= fo_ho
    0 <= fo_lo
    0 <= fo_co
    0 <= fo_res
    0 <= cracked_lo
    0 <= cracked_ho
GENERAL
    d1
    d2
END
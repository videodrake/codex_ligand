reinitialize
load cluster_results/C01_M01_S-12.45.pdb
load cluster_results/C02_M01_S-12.28.pdb
load cluster_results/C03_M01_S-11.73.pdb
load cluster_results/C04_M01_S-11.39.pdb
load cluster_results/C05_M01_S-10.37.pdb
load cluster_results/C06_M01_S-10.26.pdb
load cluster_results/C07_M01_S-10.18.pdb
load cluster_results/C08_M01_S-10.16.pdb
load cluster_results/C09_M01_S-10.12.pdb
load cluster_results/C09_M02_S-8.92.pdb
load cluster_results/C10_M01_S-9.91.pdb
load cluster_results/C11_M01_S-9.86.pdb
load cluster_results/C12_M01_S-9.81.pdb
load cluster_results/C13_M01_S-9.10.pdb
load cluster_results/C14_M01_S-9.06.pdb
load cluster_results/C15_M01_S-8.96.pdb
load cluster_results/C16_M01_S-8.91.pdb
load cluster_results/C17_M01_S-8.70.pdb
load cluster_results/C18_M01_S-8.66.pdb
load cluster_results/C19_M01_S-8.65.pdb
load cluster_results/C20_M01_S-8.42.pdb
load cluster_results/C21_M01_S-8.37.pdb
load cluster_results/C22_M01_S-8.33.pdb
load cluster_results/C23_M01_S-8.30.pdb
load cluster_results/C24_M01_S-8.22.pdb
load cluster_results/C25_M01_S-7.99.pdb
load cluster_results/C26_M01_S-7.92.pdb
load cluster_results/C27_M01_S-7.90.pdb
load cluster_results/C28_M01_S-7.85.pdb
load cluster_results/C29_M01_S-7.81.pdb
load cluster_results/C30_M01_S-7.56.pdb
load cluster_results/C31_M01_S-7.45.pdb
load cluster_results/C32_M01_S-7.26.pdb
load cluster_results/C33_M01_S-7.14.pdb
load cluster_results/C34_M01_S-7.12.pdb
load cluster_results/C35_M01_S-7.07.pdb
python
first = cmd.get_object_list()[0]
for obj in cmd.get_object_list()[1:]:
    cmd.align(f'{obj} and chain A', f'{first} and chain A')
python end
hide all
show cartoon
color gray80, chain A
group C01, C01_*
color red, C01_* and chain B
group C02, C02_*
color blue, C02_* and chain B
group C03, C03_*
color green, C03_* and chain B
group C04, C04_*
color orange, C04_* and chain B
group C05, C05_*
color cyan, C05_* and chain B
group C06, C06_*
color magenta, C06_* and chain B
group C07, C07_*
color yellow, C07_* and chain B
group C08, C08_*
color purple, C08_* and chain B
group C09, C09_*
color salmon, C09_* and chain B
group C10, C10_*
color lime, C10_* and chain B
group C11, C11_*
color slate, C11_* and chain B
group C12, C12_*
color hotpink, C12_* and chain B
group C13, C13_*
color olive, C13_* and chain B
group C14, C14_*
color teal, C14_* and chain B
group C15, C15_*
color violet, C15_* and chain B
group C16, C16_*
color brown, C16_* and chain B
group C17, C17_*
color pink, C17_* and chain B
group C18, C18_*
color aquamarine, C18_* and chain B
group C19, C19_*
color gold, C19_* and chain B
group C20, C20_*
color wheat, C20_* and chain B
group C21, C21_*
color red, C21_* and chain B
group C22, C22_*
color blue, C22_* and chain B
group C23, C23_*
color green, C23_* and chain B
group C24, C24_*
color orange, C24_* and chain B
group C25, C25_*
color cyan, C25_* and chain B
group C26, C26_*
color magenta, C26_* and chain B
group C27, C27_*
color yellow, C27_* and chain B
group C28, C28_*
color purple, C28_* and chain B
group C29, C29_*
color salmon, C29_* and chain B
group C30, C30_*
color lime, C30_* and chain B
group C31, C31_*
color slate, C31_* and chain B
group C32, C32_*
color hotpink, C32_* and chain B
group C33, C33_*
color olive, C33_* and chain B
group C34, C34_*
color teal, C34_* and chain B
group C35, C35_*
color violet, C35_* and chain B

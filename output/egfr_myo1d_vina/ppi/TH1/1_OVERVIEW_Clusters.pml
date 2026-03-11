reinitialize
load cluster_results/C01_M01_S-15.47.pdb
load cluster_results/C02_M01_S-14.64.pdb
load cluster_results/C03_M01_S-14.56.pdb
load cluster_results/C03_M02_S-9.31.pdb
load cluster_results/C04_M01_S-14.07.pdb
load cluster_results/C05_M01_S-13.99.pdb
load cluster_results/C05_M02_S-6.09.pdb
load cluster_results/C06_M01_S-12.80.pdb
load cluster_results/C07_M01_S-12.33.pdb
load cluster_results/C07_M02_S-7.07.pdb
load cluster_results/C08_M01_S-12.26.pdb
load cluster_results/C08_M02_S-6.16.pdb
load cluster_results/C09_M01_S-11.85.pdb
load cluster_results/C10_M01_S-11.54.pdb
load cluster_results/C11_M01_S-11.51.pdb
load cluster_results/C12_M01_S-11.31.pdb
load cluster_results/C12_M02_S-8.89.pdb
load cluster_results/C13_M01_S-11.14.pdb
load cluster_results/C13_M02_S-7.24.pdb
load cluster_results/C14_M01_S-10.77.pdb
load cluster_results/C15_M01_S-10.51.pdb
load cluster_results/C15_M02_S-8.26.pdb
load cluster_results/C16_M01_S-10.34.pdb
load cluster_results/C17_M01_S-10.27.pdb
load cluster_results/C18_M01_S-10.19.pdb
load cluster_results/C19_M01_S-10.18.pdb
load cluster_results/C20_M01_S-9.82.pdb
load cluster_results/C20_M02_S-6.01.pdb
load cluster_results/C21_M01_S-9.58.pdb
load cluster_results/C22_M01_S-9.54.pdb
load cluster_results/C23_M01_S-9.53.pdb
load cluster_results/C24_M01_S-9.49.pdb
load cluster_results/C25_M01_S-9.44.pdb
load cluster_results/C26_M01_S-9.44.pdb
load cluster_results/C26_M02_S-8.55.pdb
load cluster_results/C27_M01_S-9.40.pdb
load cluster_results/C27_M02_S-8.27.pdb
load cluster_results/C28_M01_S-9.39.pdb
load cluster_results/C29_M01_S-9.34.pdb
load cluster_results/C30_M01_S-9.27.pdb
load cluster_results/C31_M01_S-9.22.pdb
load cluster_results/C31_M02_S-6.59.pdb
load cluster_results/C32_M01_S-9.17.pdb
load cluster_results/C33_M01_S-9.14.pdb
load cluster_results/C34_M01_S-9.10.pdb
load cluster_results/C34_M02_S-8.46.pdb
load cluster_results/C34_M03_S-6.87.pdb
load cluster_results/C35_M01_S-9.06.pdb
load cluster_results/C36_M01_S-8.67.pdb
load cluster_results/C37_M01_S-8.66.pdb
load cluster_results/C37_M02_S-8.07.pdb
load cluster_results/C38_M01_S-8.60.pdb
load cluster_results/C38_M02_S-6.92.pdb
load cluster_results/C39_M01_S-8.58.pdb
load cluster_results/C40_M01_S-8.54.pdb
load cluster_results/C41_M01_S-8.51.pdb
load cluster_results/C42_M01_S-8.51.pdb
load cluster_results/C43_M01_S-8.49.pdb
load cluster_results/C44_M01_S-8.48.pdb
load cluster_results/C45_M01_S-8.43.pdb
load cluster_results/C46_M01_S-8.39.pdb
load cluster_results/C47_M01_S-8.34.pdb
load cluster_results/C48_M01_S-8.33.pdb
load cluster_results/C49_M01_S-8.33.pdb
load cluster_results/C50_M01_S-8.25.pdb
load cluster_results/C51_M01_S-8.22.pdb
load cluster_results/C52_M01_S-8.22.pdb
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
group C36, C36_*
color brown, C36_* and chain B
group C37, C37_*
color pink, C37_* and chain B
group C38, C38_*
color aquamarine, C38_* and chain B
group C39, C39_*
color gold, C39_* and chain B
group C40, C40_*
color wheat, C40_* and chain B
group C41, C41_*
color red, C41_* and chain B
group C42, C42_*
color blue, C42_* and chain B
group C43, C43_*
color green, C43_* and chain B
group C44, C44_*
color orange, C44_* and chain B
group C45, C45_*
color cyan, C45_* and chain B
group C46, C46_*
color magenta, C46_* and chain B
group C47, C47_*
color yellow, C47_* and chain B
group C48, C48_*
color purple, C48_* and chain B
group C49, C49_*
color salmon, C49_* and chain B
group C50, C50_*
color lime, C50_* and chain B
group C51, C51_*
color slate, C51_* and chain B
group C52, C52_*
color hotpink, C52_* and chain B

reinitialize
load ../../Rank01_C01_M01_S-12.45.pdb
hide all
show cartoon
spectrum b, blue_white_red, min=-5, max=5
python
first = cmd.get_object_list()[0]
for obj in cmd.get_object_list()[1:]:
    cmd.align(f'{obj} and chain A', f'{first} and chain A')
python end
select iface_A, chain A and resi 1940+1941+1945+1980
show sticks, iface_A
color cyan, iface_A
select iface_B, chain B and resi 972+973+994+995+996
show sticks, iface_B
color yellow, iface_B
deselect

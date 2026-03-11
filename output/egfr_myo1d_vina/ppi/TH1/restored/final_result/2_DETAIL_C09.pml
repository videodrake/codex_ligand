reinitialize
load ../../Rank09_C09_M01_S-11.85.pdb
hide all
show cartoon
spectrum b, blue_white_red, min=-5, max=5
python
first = cmd.get_object_list()[0]
for obj in cmd.get_object_list()[1:]:
    cmd.align(f'{obj} and chain A', f'{first} and chain A')
python end
select iface_A, chain A and resi 1953+1954+1956
show sticks, iface_A
color cyan, iface_A
select iface_B, chain B and resi 802+803+822+823
show sticks, iface_B
color yellow, iface_B
deselect

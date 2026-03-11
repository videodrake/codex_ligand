reinitialize
load ../../Rank01_C01_M01_S-15.47.pdb
hide all
show cartoon
spectrum b, blue_white_red, min=-5, max=5
python
first = cmd.get_object_list()[0]
for obj in cmd.get_object_list()[1:]:
    cmd.align(f'{obj} and chain A', f'{first} and chain A')
python end
select iface_A, chain A and resi 701+706+756+757+760+1989+1990+1991
show sticks, iface_A
color cyan, iface_A
select iface_B, chain B and resi 820+822+937+1006
show sticks, iface_B
color yellow, iface_B
deselect

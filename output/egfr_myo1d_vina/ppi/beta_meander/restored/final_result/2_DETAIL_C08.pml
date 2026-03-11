reinitialize
load ../../Rank08_C08_M01_S-10.16.pdb
hide all
show cartoon
spectrum b, blue_white_red, min=-5, max=5
python
first = cmd.get_object_list()[0]
for obj in cmd.get_object_list()[1:]:
    cmd.align(f'{obj} and chain A', f'{first} and chain A')
python end
select iface_A, chain A and resi 991+992+993+1701+1702+1706
show sticks, iface_A
color cyan, iface_A
select iface_B, chain B and resi 968+969+970+971+979
show sticks, iface_B
color yellow, iface_B
deselect

reinitialize
load ../../Rank12_C12_M01_S-11.31.pdb
hide all
show cartoon
spectrum b, blue_white_red, min=-5, max=5
python
first = cmd.get_object_list()[0]
for obj in cmd.get_object_list()[1:]:
    cmd.align(f'{obj} and chain A', f'{first} and chain A')
python end
select iface_A, chain A and resi 1750+1752+1784+1785
show sticks, iface_A
color cyan, iface_A
select iface_B, chain B and resi 902+956+957+958
show sticks, iface_B
color yellow, iface_B
deselect

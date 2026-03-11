reinitialize
load ../../Rank02_C02_M01_S-14.64.pdb
hide all
show cartoon
spectrum b, blue_white_red, min=-5, max=5
python
first = cmd.get_object_list()[0]
for obj in cmd.get_object_list()[1:]:
    cmd.align(f'{obj} and chain A', f'{first} and chain A')
python end
select iface_A, chain A and resi 1701+1706+1708+1709+1753+1756+1757+1781+1782+1783
show sticks, iface_A
color cyan, iface_A
select iface_B, chain B and resi 812+813+815+816+911+912+913+964+966
show sticks, iface_B
color yellow, iface_B
deselect

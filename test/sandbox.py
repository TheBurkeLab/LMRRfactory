import cantera as ct
mech_obj = ct.Solution('test/data/test_mech.yaml')
print(mech_obj.input_data)
# for reaction in mech_obj.reactions():
#     print(reaction)
#     print(f"  {reaction.reaction_type}")
#     print(f"  {reaction.input_data}")
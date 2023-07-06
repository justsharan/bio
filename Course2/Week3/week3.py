with open('integer_mass_table.txt', 'r') as f:
    masses = [line.strip().split() for line in f]
    INTEGER_MASSES = { mass[0]: int(mass[1]) for mass in masses }

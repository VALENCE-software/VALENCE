# usage: python testing.py name_of_output name_of_example_to_compare_output_to
#        for instance: python testing.py be.out be
#
#        if you want to test the whole example directory (which will take a decently
#        long time), one option is to write a script like this:
#
#        #!/bin/bash
#
#        for i in $(ls /path/testing/test_cases/)
#        do
#            /path/valence < /path/testing/test_cases/$i > $i.out
#            python /path/testing.py $i.out $i
#        done
#
#        This compares a few output values. It doesn't catch everything,
#        but it's a reasonable start for acceptance testing.

from sys import argv

# standard tolerances
nuc_tol = 10**-10
energy_tol = 10**-8
opt_energy_tol = 10**-6

# helper functions
def compare_number( string, number1, number2, tolerance):
    if abs(number1 - number2) < tolerance:
        print("%s passed" % string)
        return 0
    else:
        print("%s failed by ~%25.16f" % (string, abs(number1 - number2)))
        return 1

# define class
class vsvb_output:

    def __init__(self, nuclear_repulsion,
                 guess_energy, calculation_converged,
                 total_energy,file_name ) :
        self.nuclear_repulsion = nuclear_repulsion
        self.guess_energy = guess_energy
        self.calculation_converged = calculation_converged
        self.total_energy = total_energy
        self.file_name = file_name
    def compare_output(self,test):
        print("  ### Comparing %s to example %s ###" % (self.file_name, test.file_name))
        print("")
        error = 0
        error += compare_number( "  nuclear repulsion comparison",
                 self.nuclear_repulsion, test.nuclear_repulsion, nuc_tol )
        error += compare_number( "  guess energy comparison",
                 self.guess_energy, test.guess_energy, energy_tol )
# only check total energy comparison if we expect it to converge
        if test.calculation_converged :
            if not self.calculation_converged : 
                error += 1
                print("  calculation did not converge when it should have")
            error += compare_number( "  total energy comparison",
                                     self.total_energy, test.total_energy,
                                     opt_energy_tol )
        else :
            print( "  No total energy comparison, since not a orbital optimization run")
        return error
            

# here are our standard tests
# results from running in serial with obs01 on tnt.alcf.anl.gov
list_of_tests = [ vsvb_output( 0, -23.1330343073501332,
                               False, 0, "b"),

                  vsvb_output( 0, -14.5729681271985356,
                               True, -14.5729681271984806, "be"),

                  vsvb_output( 0, -14.5729681272019285,
                               True, -14.5729681272044544, "be+ndf"),

                  vsvb_output( 0, -14.5763526670230465,
                               False, 0, "be-sv"),

                  vsvb_output( 0, -1612.8234489982967261,
                               True, -1630.3106199068311071, "cu+"),

                  vsvb_output( 0, -0.4999455685829771,
                               True, -0.4999455685829771, "h"),

                  vsvb_output( 0.6872431808311208, -1.1474697575959574,
                               True, -1.1474697575960091, "h2-dz"),

                  vsvb_output( 8.9081858332635520, -75.9781013035254063,
                               True, -75.9781274674631817, "h2o-vdz"),

                  vsvb_output( 8.9081858332635537, -75.9781013035254063,
                               True, -75.9781274066777996, "h2o-vdz-dem"),

                  vsvb_output( 8.9081858332635520, -75.9884707182847166,
                               True, -75.9886131586740277, "h2o-vdz-sc1"),

                  vsvb_output( 0.6872431808311208, -1.0638106695934189,
                               True, -1.0638106696841798, "h2-sz"),

                  vsvb_output( 0, -2.8615142272282172,
                               True, -2.8615142272282426, "he"),

                  vsvb_output( 0, -0.4999455685829772,
                               True, -0.4999455685829776, "h+ndf"),

                  vsvb_output( 0, -7.4324082125278643,
                               False,0,"li"),

                  vsvb_output( 0, -7.2891697664828667,
                               False,0,"li-"),

                  vsvb_output( 0.9644785830619010,-8.0014494601459774,
                               False,0,"lih"),

                  vsvb_output( 0.9644785830619010,-8.0001821398489863,
                               True,-8.0004605271591132,"lih-sv"),

                  vsvb_output( 42.3831302721518597,-79.1652117037620258,
                               False,0,"ethane"),

                  vsvb_output( 122.5642769690673930,-158.3213927294620476,
                               False,0,"ethane2"),

                  vsvb_output( 0,-14.5883977218581968,
                               True,-14.5883977221021848,"be-sc"),

                  vsvb_output( 30.3586353059260574,-198.8081610426414159,
                               False,0,"f2-scval-p2"),

                  vsvb_output( 23.4657784730843311,-109.0462252625622739,
                               False,0,"n2.sc4val-b.p2"),

                  vsvb_output( 0,-14.3173525131283927,
                               True,-14.3173525131288084,"be-scv3s+2sc"),

                  vsvb_output( 0,-13.9617965781859112,
                               True,-13.9617965781864086,"be3s2"),

                  vsvb_output( 0,-2.1365129715011397,
                               True,-2.1365129715011482,"he1s2s"),

                  vsvb_output( 0,-2.0680536967364538,
                               True,-2.0680536967364542,"he3s-1s3s"),

                  vsvb_output( 15.3319283821224985,-75.4791177322835267,
                               True,-75.4792511484224349,"c2-pt-vshf-p2"),

                  vsvb_output( 0,-7.4310688735521992,
                               True,-7.4325745060262483,"li_opt")

                ]

# here's a dictionary for our standard tests for easy lookup
test_dictionary = {}
for test in list_of_tests:
    test_dictionary[test.file_name] = test

# now the real work: looking at the test output file

# read in the file to test, and who to test against
if len(argv) < 3:
    print("Two arguments are needed: name of test output and name of file to compare to.")
    exit()

script_name, test_output, example_to_compare_to = argv

# grep all elements out of test file

# open test_output and fill a test vsvb_output object
test = vsvb_output( 0, 0, False, 0, test_output )

with open( test_output, 'r') as f:
    for line in f:
        if "nuclear repulsion" in line:
            test.nuclear_repulsion = float(line.split()[2])
#            print "%25.16f" % float(line.split()[2])
        if "guess energy" in line:
# if you're looking at mbs code
            if "in atomic units" in line :
                test.guess_energy = float(line.split()[5])
            else :
                test.guess_energy = float(line.split()[2])
        if "calculation converged" in line:
            test.calculation_converged = True
        if "total energy" in line:
            if "in atomic units" in line :
                test.total_energy = float(line.split()[5])
            else:
                test.total_energy = float(line.split()[2])

# compare output values
if example_to_compare_to in test_dictionary:
    errors = test.compare_output( test_dictionary[example_to_compare_to] )
else:
    print("%s is not in the testing dictionary" % example_to_compare_to)

# print result 
print("Number of errors:  %d" % errors)

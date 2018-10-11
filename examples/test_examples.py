# usage: python testing.py name_of_output name_of_example_to_compare_output_to
#        for instance: python testing.py be.out be
#
#        if you want to test the whole example directory (which will take a decently
#        long time), one option is to write a script like this:
#
#        #!/bin/bash
#
#        for i in $(ls /path/examples/)
#        do
#            /path/valence.exe < /path/examples/$i > $i.out
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
        if not (self.calculation_converged == test.calculation_converged):
            error += 1
            print("  calculation did (not) converge when it should(n't) have")
        error += compare_number( "  total energy comparison",
                 self.total_energy, test.total_energy, opt_energy_tol )
        return error
            

# here are our standard tests
# results from running in serial with obs01 on tnt.alcf.anl.gov
list_of_tests = [ vsvb_output( 0,    -14.5729681271981626, 
                               True, -14.5729681271981626, "be"),

                  vsvb_output( 0,    -14.5729681271981626,
                               True, -14.5729681271981626, "be.DBF"),

                  vsvb_output( 0,    -14.5763526671383907,
                               True, -14.5763526671383907, "be.sv"),

                  vsvb_output( 0,    -14.5883977218582270,
                               True, -14.5883977218582270, "be.2SC"),

                  vsvb_output( 0,    -14.3173525131283803,
                               True, -14.3173525131283803, "be2s3s.2SC"),

                  vsvb_output( 0,    -13.9617965781858562,
                               True, -13.9617965781858562, "be3s2"),

                  vsvb_output( 42.3777749240206134, -79.1739950294178385,
                               True, -79.1739950294178385, "c2h6"),

                  vsvb_output( 82.6758971084159811, -118.1767010509024658,
                               True, -118.1767010509024658, "c3h8"),

                  vsvb_output( 13.5333409938694693, -40.1730886556952527,
                               True, -40.1730886556952527, "ch4"),

                  vsvb_output( 0,    -1638.3464829508911862,
                               True, -1638.3464829508911862, "cu+.3d10"),

                  vsvb_output( 0,    -1638.3747554835470055,
                               True, -1638.3747554835470055, "cu+.3d94s1"),

                  vsvb_output( 0,    -99.4130881819624221,
                               True, -99.4130881819624221, "f-"),

                  vsvb_output( 0,    -1261.4826418871418809,
                               True, -1261.4826418871418809, "fe2+"),

                  vsvb_output( 0,    -1260.4285111866975058,
                               True, -1260.4285111866975058, "fe3+"),

                  vsvb_output( 0,    -0.4999455685829772,
                               True, -0.4999455685829772, "h"),

                  vsvb_output( 0.7074562155614478, -1.1471520768516472,
                               True, -1.1471520768516472, "h2.dz"),

                  vsvb_output( 0.7430177607974767, -1.1368465732971214,
                               True, -1.1368465732971214, "h2.sz"),

                  vsvb_output( 9.2527676629424676, -75.9853591758974005,
                               True, -75.9853591758974005, "h2o"),

                  vsvb_output( 9.2527676629424676, -75.9959067344043291,
                               True, -75.9959067344043291, "h2o.SC"),

                  vsvb_output( 0,    -2.8615142272282199,
                               True, -2.8615142272282199, "he"),

                  vsvb_output( 0,    -2.1365129715011388,
                               True, -2.1365129715011388, "he1s2s"),

                  vsvb_output( 0,    -2.0680536967364547,
                               True, -2.0680536967364547, "he1s3sT"),

                  vsvb_output( 0,    -7.4324082125278759,
                               True, -7.4324082125278759, "li"),

                  vsvb_output( 0,    -7.4272143140964362,
                               True, -7.4326245198905436, "li_opt"),

                  vsvb_output( 0.9682433201511886, -7.9795127137798501,
                               True, -7.9795127137798501, "lih.VSHF"),

                  vsvb_output( 0.9682433201511886, -7.9790401124852277,
                               True, -7.9790401124852277, "lih.SDVB"),

                  vsvb_output( 0.9682433201511886, -7.9959273242325057,
                               True, -7.9959273242325057, "lih.SCval"),

                  vsvb_output( 0.9682433201511886, -7.8782983082260758,
                               True, -7.8782983082260758, "lih.exstate"),

                  vsvb_output( 24.0445894035220533, -108.9439496018654978,
                               True, -108.9439496018654978, "n2.VSHF")

                  vsvb_output( 138.2471998809889158, -173.1023795379115882,
                               True, -173.1023795379115882, "nme3")

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

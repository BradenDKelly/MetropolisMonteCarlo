Implementation for rigid polyatomic Metropolis Monte Carlo with Ewald summation for electrostatics. 

While individual functions have been previously tested and passed, I have not thrown them all together into a cohesive simulation. So, there will be some kinks to work out.

- need to add in wolf summation

x (Completed) need to test on fixed configurations provides by NIST for SPC/E water https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10a-cutoff 

(need to add NIST SPC/E tests as unit tests that can be routinely called for testing).

Currently SPC/E works in the NVT ensemble (need to upload this code and upload a proper plot of an rdf, not the cliche excel default)

![](spce_rdf.png)



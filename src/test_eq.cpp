#include <stdio.h>
#include <mpfr.h>
#include "eq.hpp"

int main(void){
    primersim::Primeanneal pa;
    pa.num_cpu = 10;
    //double temp_c = 60;
    double dna_conc = 3e-15;
    double primer_conc = 250e-9;
    double mv = 100;
    double dv = 1.5;
    double dntp = 0.2;

    //pa.read_primers_individual("./codes.csv");
    //pa.assign_addresses("addresses_new.csv", 60, mv, dv, dntp);
    //pa.evaluate_addresses("addresses.csv", "sort_eval.csv", dna_conc, primer_conc, mv, dv, dntp);
    std::vector<double> temps;
    double temp = 60;
    for(int i = 0; i < 50; i++){
        temps.push_back(temp);
        /*
        temp -= 1.;
        if (temp < 40)
            temp = 40;
            */
    }

    pa.sim_pcr("sort_eval.csv", "pcr_60_test_new.csv", 0, temps.size(), temps, dna_conc, primer_conc, primer_conc, mv, dv, dntp);
    printf("%s\n%s\n", pa.addresses[0].f, pa.addresses[0].f_rc);
    /*
    for (double power = -20; power < -5.5; power += 1.){
        char filename[80];
        sprintf(filename, "sort_eval_new_250e%d.csv", (int) power);
        primer_conc = 250 *pow(10, power);
        pa.sim_pcr("sort_eval_new.csv", filename, 0, 30, dna_conc, primer_conc, primer_conc, mv, dv, dntp);
    }
    */
    return 0;
}
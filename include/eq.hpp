#include <vector>
#include <string>
#include "thal.h"
#include <mutex>

#ifndef EQ_HPP
#define EQ_HPP

#define FLOAT_PREC (1024)
/*
Initial concentrations only include indecies 0-5 inclusive
Final concentrations include all indecies
Concentration array indecies:
 0, F:  Forward primer
 1, R:  Reverse primer
 2, A:  Forward binding site
 3, B:  Reverse binding site
 4, X:  Nonspec binding site 1
 5, Y:  Nonspec bidning site 2
 6, FH: Forward primer hairpin
 7, RH: Reverse primer hairpin
 8, FF: Forward primer dimer
 9, RR: Reverse primer dimer
10, FR: Forward-Reverse primer dimer
11, FA: Forward specific binding
12, FB: Forward bound to reverse binding site
13, RB: Reverse specific binding
14, RA: Reverse bound to forward binding site
15, FX: Forward bound to nonspec 1
16, FY: Forward bound to nonspec 2
17, RX: Reverse bound to nonspec 1
18, RY: Reverse bound to nonspec 2

Equilibrium constant array indecies:
 0, K_FH:   Forward hairpin
 1, K_RH:   Reverse hairpin
 2, K_FF:   Forward dimer
 3, K_RR:   Reverse dimer
 4, K_FR:   Forward-Reverse dimer
 5, K_FA:   Forward specific
 6, K_FB:   Forward to reverse binding site
 7, K_RB:   Reverse specific
 8, K_RA:   Reverse to forward binding site
 9, K_FX:   Forward to nonspec 1
10, K_FY:   Forward to nonspec 2
11, K_RX:   Reverse to nonspec 1
12, K_RY:   Reverse to nonspec 2
*/

#define F   0
#define R   1
#define A   2
#define B   3
#define X   4
#define Y   5
#define FH  6
#define RH  7
#define FF  8
#define RR  9
#define FR  10
#define FA  11
#define FB  12
#define RB  13
#define RA  14
#define FX  15
#define FY  16
#define RX  17
#define RY  18

#define K_FH    0
#define K_RH    1
#define K_FF    2
#define K_RR    3
#define K_FR    4
#define K_FA    5
#define K_FB    6
#define K_RB    7
#define K_RA    8
#define K_FX    9
#define K_FY    10
#define K_RX    11
#define K_RY    12

class address_k_conc{
    public:
        mpfr_t conc_f;
        mpfr_t conc_r;
        mpfr_t conc_frc;
        mpfr_t conc_rrc;
        mpfr_t k_primer_f_addr_f;
        mpfr_t k_primer_f_addr_frc;
        mpfr_t k_primer_f_addr_r;
        mpfr_t k_primer_f_addr_rrc;
        mpfr_t k_primer_r_addr_f;
        mpfr_t k_primer_r_addr_frc;
        mpfr_t k_primer_r_addr_r;
        mpfr_t k_primer_r_addr_rrc;

        address_k_conc();
        ~address_k_conc();
};

class EQ{
    public:
        mpfr_t last_val[2];
        mpfr_t bounds[3];
        mpfr_t solutions[3];
        mpfr_t tmp[4];
        mpfr_t c[19];
        mpfr_t c0[6];
        mpfr_t k[13];
        mpfr_t spec_exp_amp;
        mpfr_t spec_lin_amp;
        mpfr_t nonspec_exp_amp;
        mpfr_t nonspec_lin_amp;
        mpfr_t best_spec_exp_amp;
        mpfr_t best_spec_lin_amp;
        mpfr_t best_nonspec_exp_amp;
        mpfr_t best_nonspec_lin_amp;
        std::vector<address_k_conc> address_k_conc_vec;

        EQ();
        ~EQ();
        void calc_ca_cb_cx_cy();
        void calc_cx();
        void calc_cy();
        void calc_cf(mpfr_t ret);
        void calc_cr(mpfr_t ret);
        void calc_bound_concs();
        void solve_cf_cr(int F_or_R);
        void solve_eq();
        void print_state(std::string out_filename, std::string s);
};

class Primeanneal{
    public:
        struct primer_info{
            char *f;
            char *rc;
            double spec_dg;
            double nonspec_dg_f;
            double nonspec_dg_rc;
            double nonspec_dg_f_f;
            double nonspec_dg_f_rc;
            double nonspec_dg_rc_f;
            double nonspec_dg_rc_rc;
        };
        struct address{
            char *f;
            char *r;
            char *f_rc;
            char *r_rc;
            double temp_c;
        };
        std::vector<Primeanneal::primer_info> primers;
        std::vector<Primeanneal::address> addresses;
        //For multithreading address evaluation
        unsigned int address_index;
        std::mutex addr_mtx;
        std::mutex outfile_mtx;
        std::mutex mpfr_mtx;
        pthread_barrier_t bar;
        unsigned int num_cpu;


    public:
        Primeanneal();
        ~Primeanneal();
        void read_primers_individual(std::string filename);
        void swap_primer_rc(Primeanneal::primer_info *p);
        void reverse_comp(const char *fwd, char *rc);
        double dg_to_eq_const(double dg, double temp_c);
        double eq_const_to_dg(double eq_const, double temp_c);
        void dg_to_eq_const_mpfr(mpfr_t ret, double dg, double temp_c);
        double eq_const_to_dg_mpfr(mpfr_t tmp, mpfr_t eq_const, double temp_c);
        void assign_addresses(const char *out_filename, double temp_c, double mv_conc, double dv_conc, double dntp_conc);
        void assign_addresses_nosort(const char *out_filename);
        //void address_access(Primeanneal::address, double temp_c, double mv_conc, double dv_conc, double dntp_conc);
        //void address_access(Primeanneal::address addr, double *f_nsrc_dh, double *f_nsrc_ds, double *r_nsrc_dh, double *r_nsrc_ds,  double *f_nsf_dh, double *f_nsf_ds, double *r_nsf_dh, double *r_nsf_ds, double temp_c, double mv_conc, double dv_conc, double dntp_conc);
        thal_results calc_dimer(char *seq1, char *seq2, thal_args a);
        thal_results calc_hairpin(char *seq1, thal_args a);
        bool primer_compare(const Primeanneal::primer_info& a, const Primeanneal::primer_info& b);
        void read_addresses(const char *filename, bool with_temp_c);
        void evaluate_addresses(const char * in_filename, const char *out_filename, double dna_conc, double primer_conc, double mv_conc, double dv_conc, double dntp_conc);
        void eval_thread(const char * out_filename, double dna_conc, double primer_conc, double mv_conc, double dv_conc, double dntp_conc);
        //void sim_pcr_thread(unsigned int tid, unsigned int num_cpu, const char *out_filename, unsigned int addr, unsigned int pcr_cycles, double mv_conc, double dv_conc, double dntp_conc);
        void sim_pcr(const char *in_filename, const char *out_filename, unsigned int addr, unsigned int pcr_cycles, double dna_conc, double primer_f_conc, double primer_r_conc, double mv_conc, double dv_conc, double dntp_conc);
};

#endif


#include <stdio.h>
#include <math.h>
#include <vector>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <mpfr.h>
#include <algorithm>
#include <thread>

#include "thal.h"
#include "eq.hpp"

namespace primersim{

    address_k_conc::address_k_conc(){
        for(int i = 0; i < 5; i++){
            for(int j = 0; j < 5; j++){
                mpfr_init2(fstrand[i][j], FLOAT_PREC);
                mpfr_init2(rstrand[i][j], FLOAT_PREC);

                mpfr_init2(fstrand_change[i][j], FLOAT_PREC);
                mpfr_init2(rstrand_change[i][j], FLOAT_PREC);
            }
        }
        mpfr_init2(last_f_conc, FLOAT_PREC);
        mpfr_init2(last_r_conc, FLOAT_PREC);
        mpfr_init2(total_f_conc, FLOAT_PREC);
        mpfr_init2(total_r_conc, FLOAT_PREC);
    }

    address_k_conc::~address_k_conc(){
        for(int i = 0; i < 5; i++){
            for(int j = 0; j < 5; j++){
                mpfr_clear(fstrand[i][j]);
                mpfr_clear(rstrand[i][j]);
                mpfr_clear(fstrand_change[i][j]);
                mpfr_clear(rstrand_change[i][j]);
            }
        }
        mpfr_clear(last_f_conc);
        mpfr_clear(last_r_conc);
        mpfr_clear(total_f_conc);
        mpfr_clear(total_r_conc);
    }

    EQ::EQ(){
        int i;
        for(i = 0; i < 2; i++)
            mpfr_init2(last_val[i], FLOAT_PREC);
        for(i = 0; i < 3; i++)
            mpfr_init2(bounds[i], FLOAT_PREC);
        for(i = 0; i < 3; i++)
            mpfr_init2(solutions[i], FLOAT_PREC);
        for(i = 0; i < 4; i++)
            mpfr_init2(tmp[i], FLOAT_PREC);
        for(i = 0; i < 19; i++)
            mpfr_init2(c[i], FLOAT_PREC);
        for(i = 0; i < 6; i++)
            mpfr_init2(c0[i], FLOAT_PREC);
        for(i = 0; i < 13; i++)
            mpfr_init2(k[i], FLOAT_PREC);
        mpfr_init2(spec_fwd_amp, FLOAT_PREC);
        mpfr_init2(spec_rev_amp, FLOAT_PREC);
        mpfr_init2(saved_frc_conc, FLOAT_PREC);
        mpfr_init2(saved_rrc_conc, FLOAT_PREC);
        mpfr_init2(nonspec_exp_amp, FLOAT_PREC);
        mpfr_init2(nonspec_lin_amp, FLOAT_PREC);
        mpfr_init2(best_spec_exp_amp, FLOAT_PREC);
        mpfr_init2(best_spec_lin_amp, FLOAT_PREC);
        mpfr_init2(best_nonspec_exp_amp, FLOAT_PREC);
        mpfr_init2(best_nonspec_lin_amp, FLOAT_PREC);
        mpfr_init2(avg_nonspec_amp, FLOAT_PREC);
        mpfr_init2(last_nonspec_frc_total, FLOAT_PREC);
        mpfr_init2(last_nonspec_rrc_total, FLOAT_PREC);
    }

    EQ::~EQ(){
        int i;
        for(i = 0; i < 2; i++)
            mpfr_clear(last_val[i]);
        for(i = 0; i < 3; i++)
            mpfr_clear(bounds[i]);
        for(i = 0; i < 3; i++)
            mpfr_clear(solutions[i]);
        for(i = 0; i < 4; i++)
            mpfr_clear(tmp[i]);
        for(i = 0; i < 19; i++)
            mpfr_clear(c[i]);
        for(i = 0; i < 6; i++)
            mpfr_clear(c0[i]);
        for(i = 0; i < 13; i++)
            mpfr_clear(k[i]);
        mpfr_clear(spec_fwd_amp);
        mpfr_clear(spec_rev_amp);
        mpfr_clear(nonspec_exp_amp);
        mpfr_clear(nonspec_lin_amp);
        mpfr_clear(best_spec_exp_amp);
        mpfr_clear(best_spec_lin_amp);
        mpfr_clear(best_nonspec_exp_amp);
        mpfr_clear(best_nonspec_lin_amp);
        mpfr_clear(avg_nonspec_amp);
        mpfr_clear(last_nonspec_frc_total);
        mpfr_clear(last_nonspec_rrc_total);
    }

    void EQ::print_state(std::string out_filename, std::string s){
        int i;
        FILE *outfile = fopen(out_filename.c_str(), "a");
        fprintf(outfile, "%s,", s.c_str());
        for(i = 0; i < 2; i++)
            mpfr_fprintf(outfile, "last_val[%d],%.9Re,", i, last_val[i]);
        for(i = 0; i < 3; i++)
            mpfr_fprintf(outfile, "bounds[%d],%.9Re,",i,bounds[i]);
        for(i = 0; i < 3; i++)
            mpfr_fprintf(outfile, "solutions[%d],%.9Re,",i,solutions[i]);
        for(i = 0; i < 4; i++)
            mpfr_fprintf(outfile, "tmp[%d],%.9Re,",i,tmp[i]);
        for(i = 0; i < 19; i++)
            mpfr_fprintf(outfile, "c[%d],%.9Re,",i,c[i]);
        for(i = 0; i < 6; i++)
            mpfr_fprintf(outfile, "c0[%d],%.9Re,",i,c0[i]);
        for(i = 0; i < 13; i++)
            mpfr_fprintf(outfile, "k[%d],%.9Re,",i,k[i]);
        mpfr_fprintf(outfile, "spec_exp,%.9Re,",spec_fwd_amp);
        mpfr_fprintf(outfile, "spec_lin,%.9Re,",spec_rev_amp);
        mpfr_fprintf(outfile, "nonspec_exp,%.9Re,",nonspec_exp_amp);
        mpfr_fprintf(outfile, "nonspec_lin%.9Re,",nonspec_lin_amp);
        mpfr_fprintf(outfile, "best_spec_exp,%.9Re,",best_spec_exp_amp);
        mpfr_fprintf(outfile, "best_spec_lin,%.9Re,",best_spec_lin_amp);
        mpfr_fprintf(outfile, "best_nonspec_exp,%.9Re,",best_nonspec_exp_amp);
        mpfr_fprintf(outfile, "best_nonspec_lin,%.9Re\n",best_nonspec_lin_amp);
        fclose(outfile);
    }

    Primeanneal::Primeanneal(){
        num_cpu = 1;
    }

    Primeanneal::~Primeanneal(){
        for(auto &p : primers){
            free(p.f);
            free(p.rc);
        }
        primers.clear();
    }

    double Primeanneal::dg_to_eq_const(double dg, double temp_c){
        double temp_k = temp_c + 273.15;
        const double ideal_gas_constant = 1.98720425864083;
        return exp(dg / (-ideal_gas_constant*temp_k));
    }

    double Primeanneal::eq_const_to_dg(double eq_const, double temp_c){
        double temp_k = temp_c + 273.15;
        const double ideal_gas_constant = 1.98720425864083;
        return log(eq_const) * (-ideal_gas_constant*temp_k);
    }

    void Primeanneal::dg_to_eq_const_mpfr(mpfr_t ret, double dg, double temp_c){
        double temp_k = temp_c + 273.15;
        mpfr_set_d(ret, dg / (-1.98720425864083 * temp_k), MPFR_RNDN);
        mpfr_exp(ret, ret, MPFR_RNDN);
    }

    
    double Primeanneal::dhds_to_eq_const(double dhds[2], double temp_c){
        double dg = dhds[dh] - (temp_c+273.15) * dhds[ds];
        return dg_to_eq_const(dg, temp_c);
    }

    double Primeanneal::eq_const_to_dg_mpfr(mpfr_t tmp, mpfr_t eq_const, double temp_c){
        double temp_k = temp_c + 273.15;
        mpfr_log(tmp, eq_const, MPFR_RNDN);
        mpfr_mul_d(tmp, tmp, -1.98720425864083 * temp_k, MPFR_RNDN);
        return mpfr_get_d(tmp, MPFR_RNDN);
    }

    //c[F], c[R], c0[:], and k[:] must be defined
    void EQ::calc_cx(){
        //c[X] = c0[X] / (1 + k[K_FX]*c[F] + k[K_RX]*c[R])
        mpfr_mul(tmp[2], k[K_FX], c[F], MPFR_RNDN);
        mpfr_mul(tmp[3], k[K_RX], c[R], MPFR_RNDN);
        mpfr_add(tmp[2], tmp[2], tmp[3], MPFR_RNDN);
        mpfr_add_d(tmp[2], tmp[2], 1.0, MPFR_RNDN);
        mpfr_div(c[X], c0[X], tmp[2], MPFR_RNDN);
    }

    //0 = c[R] + K[K_RH]*c[F] + 2*k[K_RR]*c[R]*c[R] + k[K_FR]*c[F]*c[R] + k[K_RA]*c[R]*c[A] + k[K_RB]*c[R]*c[B] + k[K_RX]*c[R]*c[X] + k[K_RY]*c[R]*c[Y] - c0[R]
    //0 = c[F] + k[K_FH]*c[F] + 2*k[K_FF]*c[F]*c[F] + k[K_FR]*c[F]*c[R] + k[K_FA]*c[F]*c[A] + k[K_FB]*c[F]*c[B] + k[K_FX]*c[F]*c[X] + k[K_FY]*c[F]*c[Y] - c0[F]


    //c[F], c[R], c0[:], and k[:] must be defined
    void EQ::calc_cf(mpfr_t ret){
        //tmp[0] = c[F] + k[K_FH]*c[F] + 2*k[K_FF]*c[F]*c[F] + k[K_FR]*c[F]*c[R] + k[K_FX]*c[F]*c[X] - c0[F]
        //tmp[0] = c[F]
        mpfr_set(tmp[0], c[F], MPFR_RNDN);
        //tmp[0] += k[K_FH]*c[F]
        mpfr_mul(tmp[2], k[K_FH], c[F], MPFR_RNDN);
        mpfr_add(tmp[0], tmp[0], tmp[2], MPFR_RNDN);
        //tmp[0] += 2*k[K_FF]*c[F]*c[F]
        mpfr_mul_d(tmp[2], k[K_FF], 2.0, MPFR_RNDN);
        mpfr_mul(tmp[2], tmp[2], c[F], MPFR_RNDN);
        mpfr_mul(tmp[2], tmp[2], c[F], MPFR_RNDN);
        mpfr_add(tmp[0], tmp[0], tmp[2], MPFR_RNDN);
        //tmp[0] += k[K_FR]*c[F]*c[R]
        mpfr_mul(tmp[2], k[K_FR], c[F], MPFR_RNDN);
        mpfr_mul(tmp[2], tmp[2], c[R], MPFR_RNDN);
        mpfr_add(tmp[0], tmp[0], tmp[2], MPFR_RNDN);
        //tmp[0] += k[K_FX]*c[F]*c[X]
        mpfr_mul(tmp[2], k[K_FX], c[F], MPFR_RNDN);
        mpfr_mul(tmp[2], tmp[2], c[X], MPFR_RNDN);
        mpfr_add(tmp[0], tmp[0], tmp[2], MPFR_RNDN);
        //ret = tmp[0] - c0[F]
        mpfr_sub(ret, tmp[0], c0[F], MPFR_RNDN);
    }

    //c[F], c[R], c0[:], and k[:] must be defined
    void EQ::calc_cr(mpfr_t ret){
        //tmp[1] = c[R] + K[K_RH]*c[F] + 2*k[K_RR]*c[R]*c[R] + k[K_FR]*c[F]*c[R] + k[K_RX]*c[R]*c[X] - c0[R]
        //tmp[1] = c[R]
        mpfr_set(tmp[1], c[R], MPFR_RNDN);
        //tmp[1] += k[K_RH]*c[R]
        mpfr_mul(tmp[2], k[K_RH], c[R], MPFR_RNDN);
        mpfr_add(tmp[1], tmp[1], tmp[2], MPFR_RNDN);
        //tmp[1] += 2*k[K_RR]*c[R]*c[R]
        mpfr_mul_d(tmp[2], k[K_RR], 2.0, MPFR_RNDN);
        mpfr_mul(tmp[2], tmp[2], c[R], MPFR_RNDN);
        mpfr_mul(tmp[2], tmp[2], c[R], MPFR_RNDN);
        mpfr_add(tmp[1], tmp[1], tmp[2], MPFR_RNDN);
        //tmp[0] += k[K_FR]*c[F]*c[R]
        mpfr_mul(tmp[2], k[K_FR], c[F], MPFR_RNDN);
        mpfr_mul(tmp[2], tmp[2], c[R], MPFR_RNDN);
        mpfr_add(tmp[0], tmp[0], tmp[2], MPFR_RNDN);
        //tmp[1] += k[K_RX]*c[R]*c[X]
        mpfr_mul(tmp[2], k[K_RX], c[R], MPFR_RNDN);
        mpfr_mul(tmp[2], tmp[2], c[X], MPFR_RNDN);
        mpfr_add(tmp[1], tmp[1], tmp[2], MPFR_RNDN);
        //ret = tmp[1] - c0[R]
        mpfr_sub(ret, tmp[1], c0[R], MPFR_RNDN);
    }

    void EQ::calc_bound_concs(){
        //c[FH] = k[K_FH] * c[F]
        mpfr_mul(c[FH], k[K_FH], c[F], MPFR_RNDN);

        //c[RH] = k[K_RH] * c[R]
        mpfr_mul(c[RH], k[K_RH], c[R], MPFR_RNDN);

        //c[FF] = k[K_FF] * c[F] * c[F]
        mpfr_mul(c[FF], k[K_FF], c[F], MPFR_RNDN);
        mpfr_mul(c[FF], c[FF], c[F], MPFR_RNDN);

        //c[RR] = k[K_RR] * c[R] * c[R]
        mpfr_mul(c[RR], k[K_RR], c[R], MPFR_RNDN);
        mpfr_mul(c[RR], c[RR], c[R], MPFR_RNDN);

        //c[FR] = k[K_FR] * c[F] * c[R]
        mpfr_mul(c[FR], k[K_FR], c[F], MPFR_RNDN);
        mpfr_mul(c[FR], c[FR], c[R], MPFR_RNDN);

        //c[FX] = k[K_FX] * c[F] * c[X]
        mpfr_mul(c[FX], k[K_FX], c[F], MPFR_RNDN);
        mpfr_mul(c[FX], c[FX], c[X], MPFR_RNDN);

        //c[RX] = k[K_RX] * c[T] * c[X]
        mpfr_mul(c[RX], k[K_RX], c[R], MPFR_RNDN);
        mpfr_mul(c[RX], c[RX], c[X], MPFR_RNDN);
    }

    /*
    bounds[0] = lower bound
    bounds[1] = midpoint
    bounds[2] = upper bound

    Solutions are the value targeting zero
    The closer to zero the more accurate the concentrations
    solutions[0] = lower bound solution
    solutions[1] = midpoint solution
    solutions[2] = upper bound solution
    */
    void EQ::solve_cf_cr(int F_or_R){
        mpfr_set_d(bounds[0], 0.0, MPFR_RNDN);
        mpfr_set(bounds[2], c0[F_or_R], MPFR_RNDN);
        mpfr_div_d(bounds[1], bounds[2], 2.0, MPFR_RNDN);

        mpfr_set(c[F_or_R], bounds[0], MPFR_RNDN);
        calc_cx();
        if(F_or_R == F)
            calc_cf(solutions[0]);
        else
            calc_cr(solutions[0]);
        mpfr_set(c[F_or_R], bounds[1], MPFR_RNDN);
        calc_cx();
        if(F_or_R == F)
            calc_cf(solutions[1]);
        else
            calc_cr(solutions[1]);
        mpfr_set(c[F_or_R], bounds[2], MPFR_RNDN);
        calc_cx();
        if(F_or_R == F)
            calc_cf(solutions[2]);
        else
            calc_cr(solutions[2]);

        //if ((upper_sol > 0.0 && mid_sol > 0.0 && lower_sol > 0.0) || (upper_sol < 0.0 && mid_sol < 0.0 && lower_sol < 0.0)){
        if((mpfr_cmp_d(solutions[2], 0.) > 0 && mpfr_cmp_d(solutions[1], 0.) > 0 && mpfr_cmp_d(solutions[0], 0.) > 0) ||
            (mpfr_cmp_d(solutions[2], 0.) < 0 && mpfr_cmp_d(solutions[1], 0.) < 0 && mpfr_cmp_d(solutions[0], 0.) < 0)){
            printf("No Solution. Needs more precision\n");
            return;
        }

        //Binary search to find the root
        for (int i = 0; i < FLOAT_PREC * 2; i++){
            //Check if root has been found
            //if (upper_sol == zero) return upper_bound;
            if(mpfr_cmp_d(solutions[2], 0.) == 0){
                mpfr_set(c[F_or_R], bounds[2], MPFR_RNDN);
                break;
            }

            //if (mid_sol == zero) return midpoint;
            if(mpfr_cmp_d(solutions[1], 0.) == 0){
                mpfr_set(c[F_or_R], bounds[1], MPFR_RNDN);
                break;
            }

            //if (lower_sol == zero) return lower_bound;
            if(mpfr_cmp_d(solutions[0], 0.) == 0){
                mpfr_set(c[F_or_R], bounds[0], MPFR_RNDN);
                break;
            }

            //Find next window
            //if ((upper_sol > zero && mid_sol < zero) || (upper_sol < zero && mid_sol > zero)){
            if((mpfr_cmp_d(solutions[2], 0.) > 0 && mpfr_cmp_d(solutions[1], 0.) < 0) || (mpfr_cmp_d(solutions[2], 0.) < 0 && mpfr_cmp_d(solutions[1], 0.) > 0)){
                //lower_bound = midpoint;
                mpfr_set(bounds[0], bounds[1], MPFR_RNDN);
                //midpoint = lower_bound / two + upper_bound / two;
                mpfr_div_d(tmp[0], bounds[0], 2., MPFR_RNDN);
                mpfr_div_d(bounds[1], bounds[2], 2., MPFR_RNDN);
                mpfr_add(bounds[1], bounds[1], tmp[0], MPFR_RNDN);
                //lower_sol = mid_sol;
                mpfr_set(solutions[0], solutions[1], MPFR_RNDN);
            } else {    //c[X] = c0[X] / (1 + k[K_FX]*c[F] + k[K_RX]*c[R])
        mpfr_mul(tmp[2], k[K_FX], c[F], MPFR_RNDN);
        mpfr_mul(tmp[3], k[K_RX], c[R], MPFR_RNDN);
        mpfr_add(tmp[2], tmp[2], tmp[3], MPFR_RNDN);
        mpfr_add_d(tmp[2], tmp[2], 1.0, MPFR_RNDN);
        mpfr_div(c[X], c0[X], tmp[2], MPFR_RNDN);
                //upper_bound = midpoint;
                mpfr_set(bounds[2], bounds[1], MPFR_RNDN);
                //midpoint = lower_bound / two + upper_bound / two;
                mpfr_div_d(tmp[0], bounds[0], 2., MPFR_RNDN);
                mpfr_div_d(bounds[1], bounds[2], 2., MPFR_RNDN);
                mpfr_add(bounds[1], bounds[1], tmp[0], MPFR_RNDN);
                //upper_sol = mid_sol;
                mpfr_set(solutions[2], solutions[1], MPFR_RNDN);
            }

            mpfr_set(c[F_or_R], bounds[1], MPFR_RNDN);
            //Check if at the limits of precision
            //if(midpoint >= upper_bound || midpoint <= lower_bound){
            if(mpfr_cmp(bounds[1], bounds[2]) >= 0 || mpfr_cmp(bounds[1], bounds[0]) <= 0){
                break;
            }

            //Calculate solution for new midpoint
            mpfr_set(c[F_or_R], bounds[1], MPFR_RNDN);
            calc_cx();
            if(F_or_R == F)
                calc_cf(solutions[1]);
            else
                calc_cr(solutions[1]);
        }
    }

    void EQ::solve_eq(){
        mpfr_div_d(c[F], c0[F], 2.0, MPFR_RNDN);
        mpfr_div_d(c[R], c0[R], 2.0, MPFR_RNDN);
        mpfr_set(last_val[F], c[F], MPFR_RNDN);
        mpfr_set(last_val[R], c[R], MPFR_RNDN);
        int i;
        for(i = 0; i < FLOAT_PREC * 2; i++){
            solve_cf_cr(F);
            //mpfr_printf("%.5Re %.5Re %.5Re\n", solutions[0], solutions[1], solutions[2]);
            solve_cf_cr(R);
            //mpfr_printf("%.5Re %.5Re %.5Re\n\n", solutions[0], solutions[1], solutions[2]);
            if(!mpfr_cmp(c[F], last_val[F]) && !mpfr_cmp(c[R], last_val[R]))
                break;
            mpfr_set(last_val[F], c[F], MPFR_RNDN);
            mpfr_set(last_val[R], c[R], MPFR_RNDN);
        }
        calc_bound_concs();
    }

    void Primeanneal::read_primers_individual(std::string filename){
        FILE *infile = fopen(filename.c_str(), "r");
        char buffer[50];
        primers.clear();
        while(fscanf(infile, "%[ACGT]%*[^\n]\n", buffer) != EOF){
            //printf("%s\n", buffer);
            int len = strlen(buffer)+1;
            if(len != 21)
                printf("%s\n", buffer);
            Primeanneal::primer_info primer_pair;
            primer_pair.f = (char *) malloc(sizeof(char) * len);
            primer_pair.rc = (char *) malloc(sizeof(char) * len);
            strcpy(primer_pair.f, buffer);
            if(strlen(primer_pair.f) != 20)
                printf("pp:%s\n", primer_pair.f);
            reverse_comp(primer_pair.f, primer_pair.rc);
            primers.push_back(primer_pair);
            if(strlen(primers[primers.size()-1].f) != 20)
                printf("%s\n", primers[primers.size()-1].f);
            if(strlen(primers[primers.size()-1].rc) != 20)
                printf("%s\n", primers[primers.size()-1].rc);
        }
        fclose(infile);
    }

    void Primeanneal::reverse_comp(const char *fwd, char *rc){
        int end = strlen(fwd)-1;
        for (int i = 0; i <= end; i++){
            if (fwd[i] == 'A')
                rc[end-i] = 'T';
            else if (fwd[i] == 'T')
                rc[end-i] = 'A';
            else if (fwd[i] == 'G')
                rc[end-i] = 'C';
            else if (fwd[i] == 'C')
                rc[end-i] = 'G';
        }
        rc[end+1] = '\0';
    }

    void Primeanneal::swap_primer_rc(Primeanneal::primer_info *p){
        char *tmp = p->f;
        p->f = p->rc;
        p->rc = tmp;
    }

    thal_results Primeanneal::calc_dimer(char *seq1, char *seq2, thal_args a){
        unsigned char *seq1_uchar = (unsigned char *)seq1;
        unsigned char *seq2_uchar = (unsigned char *)seq2;
        thal_results o;
        thal(seq1_uchar, seq2_uchar, &a, THL_FAST, &o);
        return o;
    }

    thal_results Primeanneal::calc_hairpin(char *seq1, thal_args a){
        unsigned char *seq1_uchar = (unsigned char *)seq1;
        a.type = thal_hairpin;
        thal_results o;
        thal(seq1_uchar, seq1_uchar, &a, THL_FAST, &o);
        return o;
    }

    bool Primeanneal::primer_compare(const Primeanneal::primer_info& a, const Primeanneal::primer_info& b){
        return (a.nonspec_dg_f_rc - a.spec_dg > b.nonspec_dg_f_rc - b.spec_dg);
    }

    void Primeanneal::assign_addresses(const char *out_filename, double temp_c, double mv_conc, double dv_conc, double dntp_conc){
        thal_args a;
        set_thal_default_args(&a);
        a.mv = mv_conc;
        a.dv = dv_conc;
        a.dntp = dntp_conc;
        a.temp = temp_c + 273.15;
        int half_num_nonspec = primers.size() - 1;
        bool swap = true;
        int iterations = 0;
        while(swap){
            swap = false;
            int num_swaps = 0;
            int count = 0;
            for (auto &p : primers){
                printf("%d\n", count++);
                if(strlen(p.f) != 20)
                    printf("%s\n", p.f);
                if(strlen(p.rc) != 20)
                    printf("%s\n", p.rc);
                p.spec_dg = calc_dimer(p.f, p.rc, a).dg;
                double f_nonspec_k = 0.0;
                double rc_nonspec_k = 0.0;
                double f_f_nonspec_k = 0.;
                double f_rc_nonspec_k = 0.;
                double rc_f_nonspec_k = 0.;
                double rc_rc_nonspec_k = 0.;
                for(auto &p2 : primers){
                    if(!strcmp(p.f, p2.f))
                        continue;
                    //Fwd-Fwd
                    double nonspec_dg;// = calc_dimer(p.f, p2.f, a).dg;
                    //f_nonspec_k += dg_to_eq_const(nonspec_dg, temp_c) / num_nonspec;
                    //f_f_nonspec_k += dg_to_eq_const(nonspec_dg, temp_c) / half_num_nonspec;
                    //Fwd-RC
                    nonspec_dg = calc_dimer(p.f, p2.rc, a).dg;
                    //f_nonspec_k += dg_to_eq_const(nonspec_dg, temp_c) / num_nonspec;
                    f_rc_nonspec_k += dg_to_eq_const(nonspec_dg, temp_c) / half_num_nonspec;
                    //RC-Fwd
                    //nonspec_dg = calc_dimer(p.rc, p2.f, a).dg;
                    //rc_nonspec_k += dg_to_eq_const(nonspec_dg, temp_c) / num_nonspec;
                    //rc_f_nonspec_k += dg_to_eq_const(nonspec_dg, temp_c) / half_num_nonspec;
                    //RC-RC
                    nonspec_dg = calc_dimer(p.rc, p2.rc, a).dg;
                    //rc_nonspec_k += dg_to_eq_const(nonspec_dg, temp_c) / num_nonspec;
                    rc_rc_nonspec_k += dg_to_eq_const(nonspec_dg, temp_c) / half_num_nonspec;
                }
                p.nonspec_dg_f = eq_const_to_dg(f_nonspec_k, temp_c);
                p.nonspec_dg_rc = eq_const_to_dg(rc_nonspec_k, temp_c);
                p.nonspec_dg_f_f = eq_const_to_dg(f_f_nonspec_k, temp_c);
                p.nonspec_dg_f_rc = eq_const_to_dg(f_rc_nonspec_k, temp_c);
                p.nonspec_dg_rc_f = eq_const_to_dg(rc_f_nonspec_k, temp_c);
                p.nonspec_dg_rc_rc = eq_const_to_dg(rc_rc_nonspec_k, temp_c);
                if(p.nonspec_dg_f_rc - p.spec_dg < p.nonspec_dg_rc_rc - p.spec_dg){
                    swap_primer_rc(&p);
                    swap = true;
                    num_swaps++;
                }
            }
            iterations++;
            printf("%d: %d\n", iterations, num_swaps);
        }
        std::sort(primers.begin(), primers.end(), std::bind(&Primeanneal::primer_compare, this, std::placeholders::_1, std::placeholders::_2));
        FILE *outfile = fopen("addresses.csv", "w+");
        int num_addresses = primers.size() / 2;
        for (int i = 0; i < num_addresses; i++){
            fprintf(outfile, "%s,%s,%lf,%lf,%lf,%lf,%lf,%lf\n", primers[i].f, primers[i+num_addresses].f,primers[i].spec_dg,
                                primers[i+num_addresses].spec_dg, primers[i].nonspec_dg_f_rc, primers[i+num_addresses].nonspec_dg_f_rc,
                                primers[i].nonspec_dg_f_rc - primers[i].spec_dg,
                                primers[i+num_addresses].nonspec_dg_f_rc - primers[i+num_addresses].spec_dg);
        }
        fclose(outfile);
    }

    void Primeanneal::assign_addresses_nosort(const char *out_filename){
        FILE *outfile = fopen(out_filename, "w+");
        int num_addresses = primers.size() / 2;
        for (int i = 0; i < num_addresses; i++){
            fprintf(outfile, "%s,%s,%lf,%lf,%lf,%lf,%lf,%lf\n", primers[i].f, primers[i+num_addresses].f,primers[i].spec_dg,
                                primers[i+num_addresses].spec_dg, primers[i].nonspec_dg_f_rc, primers[i+num_addresses].nonspec_dg_f_rc,
                                primers[i].nonspec_dg_f_rc - primers[i].spec_dg,
                                primers[i+num_addresses].nonspec_dg_f_rc - primers[i+num_addresses].spec_dg);
        }
        fclose(outfile);
    }

    void Primeanneal::read_addresses(const char *filename, bool with_temp_c = false){
        FILE *infile = fopen(filename, "r");
        address a;
        char tmp_f[100];
        char tmp_r[100];
        while(true){
            if(with_temp_c){
                if(fscanf(infile, "%[ACGT],%[ACGT],%lf,%*[^\n]\n", tmp_f, tmp_r, &(a.temp_c)) == EOF)
                    break;
            } else {
                if(fscanf(infile, "%[ACGT],%[ACGT]%*[^\n]\n", tmp_f, tmp_r) == EOF)
                    break;
            }
            a.f = (char *)malloc(strlen(tmp_f)+1);
            a.r = (char *)malloc(strlen(tmp_f)+1);
            a.f_rc = (char *)malloc(strlen(tmp_f)+1);
            a.r_rc = (char *)malloc(strlen(tmp_f)+1);
            strcpy(a.f, tmp_f);
            strcpy(a.r, tmp_r);
            reverse_comp(a.f, a.f_rc);
            reverse_comp(a.r, a.r_rc);
            addresses.push_back(a);
        }
    }

    void Primeanneal::eval_thread(const char * out_filename, double dna_conc, double primer_conc, double mv_conc, double dv_conc, double dntp_conc){
        EQ eq;
        mpfr_set_d(eq.c0[F], primer_conc, MPFR_RNDN);
        mpfr_set_d(eq.c0[R], primer_conc, MPFR_RNDN);
        //mpfr_set_d(eq.c0[A], dna_conc / addresses.size(), MPFR_RNDN);
        //mpfr_set_d(eq.c0[B], dna_conc / addresses.size(), MPFR_RNDN);
        thal_results o;
        thal_args ta;
        set_thal_default_args(&ta);
        ta.mv = mv_conc;
        ta.dv = dv_conc;
        ta.dntp = dntp_conc;
        double f_hp_dh;
        double f_hp_ds;
        double r_hp_dh;
        double r_hp_ds;
        struct dh_ds_info{
            double f_f_dh;
            double f_f_ds;
            double f_r_dh;
            double f_r_ds;
            double f_frc_dh;
            double f_frc_ds;
            double f_rrc_dh;
            double f_rrc_ds;
            double r_f_dh;
            double r_f_ds;
            double r_r_dh;
            double r_r_ds;
            double r_frc_dh;
            double r_frc_ds;
            double r_rrc_dh;
            double r_rrc_ds;
        };

        std::vector<struct dh_ds_info> dh_ds;
        dh_ds.resize(addresses.size());
        while(true){
            unsigned int i = 0;
            addr_mtx.lock();
            i = address_index;
            address_index++;
            addr_mtx.unlock();
            if (i >= addresses.size()){
                break;
            }
            
            double best_temp = 0.0;
            bool first = true;
            o = calc_hairpin(addresses[i].f, ta);
            f_hp_dh = o.dh;
            f_hp_ds = o.ds;
            o = calc_hairpin(addresses[i].r, ta);
            r_hp_dh = o.dh;
            r_hp_ds = o.ds;
            for (unsigned int j = 0; j < addresses.size(); j++){
                struct dh_ds_info s;
                o = calc_dimer(addresses[i].f, addresses[j].f, ta);
                s.f_f_dh = o.dh;
                s.f_f_ds = o.ds;
                o = calc_dimer(addresses[i].f, addresses[j].r, ta);
                s.f_r_dh = o.dh;
                s.f_r_ds = o.ds;
                o = calc_dimer(addresses[i].f, addresses[j].f_rc, ta);
                s.f_frc_dh = o.dh;
                s.f_frc_ds = o.ds;
                o = calc_dimer(addresses[i].f, addresses[j].r_rc, ta);
                s.f_rrc_dh = o.dh;
                s.f_rrc_ds = o.ds;
                o = calc_dimer(addresses[i].r, addresses[j].f, ta);
                s.r_f_dh = o.dh;
                s.r_f_ds = o.ds;
                o = calc_dimer(addresses[i].r, addresses[j].r, ta);
                s.r_r_dh = o.dh;
                s.r_r_ds = o.ds;
                o = calc_dimer(addresses[i].r, addresses[j].f_rc, ta);
                s.r_frc_dh = o.dh;
                s.r_frc_ds = o.ds;
                o = calc_dimer(addresses[i].r, addresses[j].r_rc, ta);
                s.r_rrc_dh = o.dh;
                s.r_rrc_ds = o.ds;
                dh_ds[j] = s;
            }

            mpfr_set_d(eq.best_spec_exp_amp, 0, MPFR_RNDN);
            mpfr_set_d(eq.best_spec_lin_amp, 0, MPFR_RNDN);
            mpfr_set_d(eq.best_nonspec_exp_amp, 0, MPFR_RNDN);
            mpfr_set_d(eq.best_nonspec_lin_amp, 0, MPFR_RNDN);

            for(double temp_c = 40.; temp_c < 80.01; temp_c += 1.){
                mpfr_set_d(eq.c0[X], dna_conc * ((4.0 * (addresses.size() - 1.0)) + 2.0), MPFR_RNDN); //4 nonspecific binding sites per double strand DNA, plus two 3' binding sites on target DNA

                mpfr_set_d(eq.spec_fwd_amp, 0, MPFR_RNDN);
                mpfr_set_d(eq.spec_rev_amp, 0, MPFR_RNDN);
                mpfr_set_d(eq.nonspec_exp_amp, 0, MPFR_RNDN);
                mpfr_set_d(eq.nonspec_lin_amp, 0, MPFR_RNDN);

                /*
                mpfr_set_d(eq.k[K_FA], dg_to_eq_const(dh_ds[i].f_frc_dh - (dh_ds[i].f_frc_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                mpfr_set_d(eq.k[K_FB], dg_to_eq_const(dh_ds[i].f_rrc_dh - (dh_ds[i].f_rrc_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                mpfr_set_d(eq.k[K_RA], dg_to_eq_const(dh_ds[i].r_frc_dh - (dh_ds[i].r_frc_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                mpfr_set_d(eq.k[K_RB], dg_to_eq_const(dh_ds[i].r_rrc_dh - (dh_ds[i].r_rrc_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                */

                mpfr_set_d(eq.k[K_FF], dg_to_eq_const(dh_ds[i].f_f_dh - (dh_ds[i].f_f_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                mpfr_set_d(eq.k[K_RR], dg_to_eq_const(dh_ds[i].r_r_dh - (dh_ds[i].r_r_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                mpfr_set_d(eq.k[K_FR], dg_to_eq_const(dh_ds[i].f_r_dh - (dh_ds[i].f_r_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);

                mpfr_set_d(eq.k[K_FH], dg_to_eq_const(f_hp_dh - (f_hp_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                mpfr_set_d(eq.k[K_RH], dg_to_eq_const(r_hp_dh - (r_hp_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);

                /*
                mpfr_set_d(eq.k[K_FY], 0., MPFR_RNDN);
                mpfr_set_d(eq.k[K_RY], 0., MPFR_RNDN);
                */
                mpfr_set_d(eq.k[K_FX], 0., MPFR_RNDN);
                mpfr_set_d(eq.k[K_RX], 0., MPFR_RNDN);
                int count = 0;
                for(unsigned int j = 0; j < addresses.size(); j++){
                    mpfr_add_d(eq.k[K_FX], eq.k[K_FX], dg_to_eq_const(dh_ds[j].f_f_dh - (dh_ds[j].f_f_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                    mpfr_add_d(eq.k[K_RX], eq.k[K_RX], dg_to_eq_const(dh_ds[j].r_f_dh - (dh_ds[j].r_f_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                    count++;

                    mpfr_add_d(eq.k[K_FX], eq.k[K_FX], dg_to_eq_const(dh_ds[j].f_frc_dh - (dh_ds[j].f_frc_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                    mpfr_add_d(eq.k[K_RX], eq.k[K_RX], dg_to_eq_const(dh_ds[j].r_frc_dh - (dh_ds[j].r_frc_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                    count++;

                    mpfr_add_d(eq.k[K_FX], eq.k[K_FX], dg_to_eq_const(dh_ds[j].f_r_dh - (dh_ds[j].f_r_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                    mpfr_add_d(eq.k[K_RX], eq.k[K_RX], dg_to_eq_const(dh_ds[j].r_r_dh - (dh_ds[j].r_r_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                    count++;

                    mpfr_add_d(eq.k[K_FX], eq.k[K_FX], dg_to_eq_const(dh_ds[j].f_rrc_dh - (dh_ds[j].f_rrc_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                    mpfr_add_d(eq.k[K_RX], eq.k[K_RX], dg_to_eq_const(dh_ds[j].r_rrc_dh - (dh_ds[j].r_rrc_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                    count++;

                }
                mpfr_div_d(eq.k[K_FX], eq.k[K_FX], count, MPFR_RNDN);
                mpfr_div_d(eq.k[K_RX], eq.k[K_RX], count, MPFR_RNDN);
                eq.solve_eq();
                //c[F] and c[R] are now solved, all nonspec concs can now be solved individually

                //TODO: Fix this
                /*
                mpfr_add(eq.tmp[0], eq.c[FA], eq.c[RA], MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.c[FB], eq.c[RB], MPFR_RNDN);
                mpfr_min(eq.spec_fwd_amp, eq.tmp[0], eq.tmp[1], MPFR_RNDN);
                mpfr_max(eq.tmp[0], eq.tmp[0], eq.tmp[1], MPFR_RNDN);
                mpfr_sub(eq.spec_rev_amp, eq.tmp[0], eq.spec_fwd_amp, MPFR_RNDN);
                */

                //Set c[X] and c[Y] to concs of individual strands
                mpfr_set_d(eq.c0[X], dna_conc / addresses.size(), MPFR_RNDN);
                //mpfr_set_d(eq.c0[Y], dna_conc / addresses.size(), MPFR_RNDN);
                for (unsigned int j = 0; j < addresses.size(); j++){
                    if (i == j)
                        continue;
                    mpfr_set_d(eq.k[K_FX], dg_to_eq_const(dh_ds[j].f_frc_dh - (dh_ds[j].f_frc_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                    mpfr_set_d(eq.k[K_RX], dg_to_eq_const(dh_ds[j].r_frc_dh - (dh_ds[j].r_frc_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                    //mpfr_set_d(eq.k[K_FY], dg_to_eq_const(dh_ds[j].f_rrc_dh - (dh_ds[j].f_rrc_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                    //mpfr_set_d(eq.k[K_RY], dg_to_eq_const(dh_ds[j].r_rrc_dh - (dh_ds[j].r_rrc_ds * (temp_c + 273.15)), temp_c), MPFR_RNDN);
                    eq.calc_cx();
                    eq.calc_bound_concs();

                    //nonspec_exp_amp += min(fwd bound, rev bound)
                    mpfr_add(eq.tmp[0], eq.c[FX], eq.c[RX], MPFR_RNDN);
                    //mpfr_add(eq.tmp[1], eq.c[FY], eq.c[RY], MPFR_RNDN);
                    mpfr_min(eq.tmp[2], eq.tmp[0], eq.tmp[1], MPFR_RNDN);
                    mpfr_add(eq.nonspec_exp_amp, eq.nonspec_exp_amp, eq.tmp[2], MPFR_RNDN);

                    //nonspec_lin_amp += max(fwd bound, rev bound) - min(fwd bound, rev bound)
                    mpfr_max(eq.tmp[0], eq.tmp[0], eq.tmp[1], MPFR_RNDN);
                    mpfr_sub(eq.tmp[2], eq.tmp[0], eq.tmp[2], MPFR_RNDN);
                    mpfr_add(eq.nonspec_lin_amp, eq.nonspec_lin_amp, eq.tmp[2], MPFR_RNDN);
                }
                mpfr_div(eq.tmp[0], eq.spec_fwd_amp, eq.nonspec_exp_amp, MPFR_RNDN);
                mpfr_div(eq.tmp[1], eq.best_spec_exp_amp, eq.best_nonspec_exp_amp, MPFR_RNDN);
                if(first || (mpfr_cmp(eq.tmp[0], eq.tmp[1]) > 0)){
                    first = false;
                    best_temp = temp_c;
                    mpfr_set(eq.best_spec_exp_amp, eq.spec_fwd_amp, MPFR_RNDN);
                    mpfr_set(eq.best_spec_lin_amp, eq.spec_rev_amp, MPFR_RNDN);
                    mpfr_set(eq.best_nonspec_exp_amp, eq.nonspec_exp_amp, MPFR_RNDN);
                    mpfr_set(eq.best_nonspec_lin_amp, eq.nonspec_lin_amp, MPFR_RNDN);
                }
            }
            mpfr_div(eq.tmp[0], eq.best_spec_exp_amp, eq.best_nonspec_exp_amp, MPFR_RNDN);
            mpfr_div(eq.tmp[1], eq.best_spec_lin_amp, eq.best_nonspec_lin_amp, MPFR_RNDN);
            outfile_mtx.lock();
            FILE *outfile = fopen(out_filename, "a");
            mpfr_fprintf(outfile, "%s,%s,%lf,%u,%.9Re,%.9Re,%.9Re,%.9Re,%.9Re,%.9Re\n",addresses[i].f, addresses[i].r, best_temp, i, eq.tmp[0], eq.tmp[1], eq.best_spec_exp_amp, eq.best_nonspec_exp_amp, eq.best_spec_lin_amp, eq.best_nonspec_lin_amp);
            fclose(outfile);
            outfile_mtx.unlock();
        }
    }

    void Primeanneal::evaluate_addresses(const char *in_filename, const char * out_filename, double dna_conc, double primer_conc, double mv_conc, double dv_conc, double dntp_conc){
        if(addresses.size() == 0)
            read_addresses(in_filename, dna_conc);
        //Clear the outfile
        FILE *outfile = fopen(out_filename, "w");
        fclose(outfile);
        std::vector<std::thread> threads;
        address_index = 0;
        for(unsigned int i = 0; i < num_cpu; i++){
            threads.push_back(std::thread(&Primeanneal::eval_thread, this, out_filename, dna_conc, primer_conc, mv_conc, dv_conc, dntp_conc));
        }

        for (auto &t : threads)
            t.join();
    }

    void Primeanneal::sim_pcr(const char * in_filename, const char *out_filename, unsigned int addr, unsigned int pcr_cycles, const std::vector<double> &temp_c_profile, double dna_conc, double primer_f_conc, double primer_r_conc, double mv_conc, double dv_conc, double dntp_conc){
        addresses.clear();
        read_addresses(in_filename, true);
        EQ eq;
        thal_args ta;
        thal_results o;
        set_thal_default_args(&ta);
        ta.mv = mv_conc;
        ta.dv = dv_conc;
        ta.dntp = dntp_conc;
        ta.temp = 273.15 + addresses[addr].temp_c;
        double temp_c = addresses[addr].temp_c;
        eq.address_k_conc_vec.resize(addresses.size());
        mpfr_set_d(eq.last_nonspec_frc_total, 0., MPFR_RNDN);
        mpfr_set_d(eq.last_nonspec_rrc_total, 0., MPFR_RNDN);

        for(unsigned int i = 0; i < addresses.size(); i++){
            mpfr_set_d(eq.address_k_conc_vec[i].fstrand[0][0], dna_conc / addresses.size(), MPFR_RNDN);
            mpfr_set_d(eq.address_k_conc_vec[i].rstrand[0][0], dna_conc / addresses.size(), MPFR_RNDN);
            for(int ii = 0; ii < 5; ii++){
                for(int jj = 0; jj < 5; jj++){
                    if (!jj && !ii)
                        continue;
                    mpfr_set_d(eq.address_k_conc_vec[i].fstrand[ii][jj], 0.0, MPFR_RNDN);
                    mpfr_set_d(eq.address_k_conc_vec[i].rstrand[ii][jj], 0.0, MPFR_RNDN);
                }
            }

            o = calc_dimer(addresses[addr].f, addresses[i].f, ta);
            eq.address_k_conc_vec[i].dhds_primer_f_addr_f[dh] = o.dh;
            eq.address_k_conc_vec[i].dhds_primer_f_addr_f[ds] = o.ds;
            //mpfr_set_d(eq.address_k_conc_vec[i].k_primer_f_addr_f, dg_to_eq_const(o.dg, temp_c), MPFR_RNDN);
            o = calc_dimer(addresses[addr].f, addresses[i].f_rc, ta);
            eq.address_k_conc_vec[i].dhds_primer_f_addr_frc[dh] = o.dh;
            eq.address_k_conc_vec[i].dhds_primer_f_addr_frc[ds] = o.ds;
            //mpfr_set_d(eq.address_k_conc_vec[i].k_primer_f_addr_frc, dg_to_eq_const(o.dg, temp_c), MPFR_RNDN);
            o = calc_dimer(addresses[addr].f, addresses[i].r, ta);
            eq.address_k_conc_vec[i].dhds_primer_f_addr_r[dh] = o.dh;
            eq.address_k_conc_vec[i].dhds_primer_f_addr_r[ds] = o.ds;
            //mpfr_set_d(eq.address_k_conc_vec[i].k_primer_f_addr_r, dg_to_eq_const(o.dg, temp_c), MPFR_RNDN);
            o = calc_dimer(addresses[addr].f, addresses[i].r_rc, ta);
            eq.address_k_conc_vec[i].dhds_primer_f_addr_rrc[dh] = o.dh;
            eq.address_k_conc_vec[i].dhds_primer_f_addr_rrc[ds] = o.ds;
            //mpfr_set_d(eq.address_k_conc_vec[i].k_primer_f_addr_rrc, dg_to_eq_const(o.dg, temp_c), MPFR_RNDN);
            o = calc_dimer(addresses[addr].r, addresses[i].f, ta);
            eq.address_k_conc_vec[i].dhds_primer_r_addr_f[dh] = o.dh;
            eq.address_k_conc_vec[i].dhds_primer_r_addr_f[ds] = o.ds;
            //mpfr_set_d(eq.address_k_conc_vec[i].k_primer_r_addr_f, dg_to_eq_const(o.dg, temp_c), MPFR_RNDN);
            o = calc_dimer(addresses[addr].r, addresses[i].f_rc, ta);
            eq.address_k_conc_vec[i].dhds_primer_r_addr_frc[dh] = o.dh;
            eq.address_k_conc_vec[i].dhds_primer_r_addr_frc[ds] = o.ds;
            //mpfr_set_d(eq.address_k_conc_vec[i].k_primer_r_addr_frc, dg_to_eq_const(o.dg, temp_c), MPFR_RNDN);
            o = calc_dimer(addresses[addr].r, addresses[i].r, ta);
            eq.address_k_conc_vec[i].dhds_primer_r_addr_r[dh] = o.dh;
            eq.address_k_conc_vec[i].dhds_primer_r_addr_r[ds] = o.ds;
            //mpfr_set_d(eq.address_k_conc_vec[i].k_primer_r_addr_r, dg_to_eq_const(o.dg, temp_c), MPFR_RNDN);
            o = calc_dimer(addresses[addr].r, addresses[i].r_rc, ta);
            eq.address_k_conc_vec[i].dhds_primer_r_addr_rrc[dh] = o.dh;
            eq.address_k_conc_vec[i].dhds_primer_r_addr_rrc[ds] = o.ds;
            //mpfr_set_d(eq.address_k_conc_vec[i].k_primer_r_addr_rrc, dg_to_eq_const(o.dg, temp_c), MPFR_RNDN);
        }

        double dhds_f_hairpin[2];
        o = calc_hairpin(addresses[addr].f, ta);
        dhds_f_hairpin[dh] = o.dh;
        dhds_f_hairpin[ds] = o.ds;

        double dhds_r_hairpin[2];
        o = calc_hairpin(addresses[addr].r, ta);
        dhds_r_hairpin[dh] = o.dh;
        dhds_r_hairpin[ds] = o.ds;

        double eq_const;
        double dg;

        mpfr_set_d(eq.c0[F], primer_f_conc, MPFR_RNDN);
        mpfr_set_d(eq.c0[R], primer_r_conc, MPFR_RNDN);

        FILE *outfile = fopen(out_filename, "w");
        //cycle,f_primer,r_primer,spec_f,spec_r,nonspec_f,nonspec_r
        mpfr_fprintf(outfile, "%02u,%lf,0.000000000,0.000000000,%.9Re,%.9Re,%.9Re,%.9Re", 0, temp_c_profile[0], eq.c0[F], eq.c0[R], eq.address_k_conc_vec[addr].fstrand[addr_end][addr_end], eq.address_k_conc_vec[addr].rstrand[addr_end][addr_end]);
        mpfr_set_d(eq.tmp[0], 0., MPFR_RNDN);
        mpfr_set_d(eq.tmp[1], 0., MPFR_RNDN);
        for(unsigned int i = 0; i < addresses.size(); i++){
            if (i == addr)
                continue;
            mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[addr_end][addr_end], MPFR_RNDN);
            mpfr_add(eq.tmp[1], eq.tmp[1], eq.address_k_conc_vec[i].rstrand[addr_end][addr_end], MPFR_RNDN);
        }
        mpfr_fprintf(outfile, ",%.9Re,0.000000000,%.9Re,0.000000000\n", eq.tmp[0], eq.tmp[1]);
        fclose(outfile);
        for(unsigned int cycle = 1; cycle <= pcr_cycles; cycle++){
            //calc c[F] and c[R]
            //eq.tmp[0] will hold the sum of the nonspecific concentrations
            //eq.tmp[1] will hold the sum of (nonspec_conc * nonspec_k_f)
            //eq.tmp[2] will hold the sum of (nonspec_conc * nonspec_k_r)
            //eq.tmp[1] / eq.tmp[0] is average nonspec k for forward primer
            //eq.tmp[2] / eq.tmp[0] is average nonspec k for reverse primer
            //eq/tmp[3] will hold partial results

            //Set the eq constant of all specific bindings, primer dimers, and primer hairpins
            eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
            mpfr_set_d(eq.k[K_FF], eq_const, MPFR_RNDN);

            eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
            mpfr_set_d(eq.k[K_FR], eq_const, MPFR_RNDN);

            eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
            mpfr_set_d(eq.k[K_RR], eq_const, MPFR_RNDN);

            eq_const = dhds_to_eq_const(dhds_f_hairpin, temp_c_profile[cycle-1]);
            mpfr_set_d(eq.k[K_FH], eq_const, MPFR_RNDN);

            eq_const = dhds_to_eq_const(dhds_r_hairpin, temp_c_profile[cycle-1]);
            mpfr_set_d(eq.k[K_RH], dg_to_eq_const(o.dg, temp_c), MPFR_RNDN);

            mpfr_set_d(eq.tmp[0], 0., MPFR_RNDN);//Total concentration
            mpfr_set_d(eq.tmp[1], 0., MPFR_RNDN);//sum(f primer eq const * partial concentration)
            mpfr_set_d(eq.tmp[2], 0., MPFR_RNDN);//sum(r primer eq const * partial concentration)

            for(unsigned int i = 0; i < addresses.size(); i++){
                //Reverse strand bindings:

                // 5' R -> R Strand -> FRC 3' [0][0]
                //5' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[addr_end][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[addr_end][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                // 5' R -> R Strand -> Primer F 3' [0][1]
                //5' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[addr_end][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[addr_end][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                // 5' R -> R Strand -> Primer FRC 3' [0][2]
                //5' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[addr_end][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[addr_end][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                // 5' R -> R Strand -> Primer R 3' [0][3]
                //5' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[addr_end][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[addr_end][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                // 5' R -> R Strand -> Primer RRC 3' [0][4]
                //5' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[addr_end][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[addr_end][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[addr_end][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer F -> R Strand -> FRC 3' [1][0]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerf][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerf][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                // 5' Primer F -> R Strand -> Primer F 3' [1][1]
                //5' and 3' Binding (same effect so no need for separate calculations)
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerf][primerf], MPFR_RNDN);
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerf][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer F -> R Strand -> Primer FRC 3' [1][2]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerf][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerf][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer F -> R Strand -> Primer R 3' [1][3]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerf][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerf][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer F -> R Strand -> Primer RRC 3' [1][4]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerf][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerf][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerf][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer FRC -> R Strand -> FRC 3' [2][0]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerfrc][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerfrc][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer FRC -> R Strand -> Primer F 3' [2][1]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerfrc][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerfrc][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer FRC -> R Strand -> Primer FRC 3' [2][2]
                //5' and 3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerfrc][primerfrc], MPFR_RNDN);
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerfrc][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer FRC -> R Strand -> Primer R 3' [2][3]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerfrc][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerfrc][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer FRC -> R Strand -> Primer RRC 3' [2][4]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerfrc][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerfrc][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerfrc][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer R -> R Strand -> FRC 3' [3][0]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerr][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerr][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer R -> R Strand -> Primer F 3' [3][1]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerr][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerr][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer R -> R Strand -> Primer FRC 3' [3][2]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerr][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerr][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer R -> R Strand -> Primer R 3' [3][3]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerr][primerr], MPFR_RNDN);
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerr][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer R -> R Strand -> Primer RRC 3' [3][4]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerr][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerr][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerr][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer RRC -> R Strand -> FRC 3' [4][0]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerrrc][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerrrc][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer RRC -> R Strand -> Primer F 3' [4][1]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerrrc][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerrrc][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer RRC -> R Strand -> Primer FRC 3' [4][2]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerrrc][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerrrc][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer RRC -> R Strand -> Primer R 3' [4][3]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerrrc][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerrrc][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer RRC -> R Strand -> Primer RRC 3' [4][4]
                //5' and 3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerrrc][primerrrc], MPFR_RNDN);
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].rstrand[primerrrc][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].rstrand[primerrrc][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //Forward Strand bindings

                // 5' R -> F Strand -> FRC 3' [0][0]
                //5' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[addr_end][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[addr_end][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                // 5' R -> F Strand -> Primer F 3' [0][1]
                //5' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[addr_end][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[addr_end][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                // 5' R -> F Strand -> Primer FRC 3' [0][2]
                //5' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[addr_end][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[addr_end][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                // 5' R -> F Strand -> Primer R 3' [0][3]
                //5' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[addr_end][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[addr_end][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                // 5' R -> F Strand -> Primer RRC 3' [0][4]
                //5' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[addr_end][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[addr_end][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[addr_end][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer F -> F Strand -> FRC 3' [1][0]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerf][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerf][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                // 5' Primer F -> F Strand -> Primer F 3' [1][1]
                //5' and 3' Binding (same effect so no need for separate calculations)
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerf][primerf], MPFR_RNDN);
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerf][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer F -> F Strand -> Primer FRC 3' [1][2]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerf][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerf][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer F -> F Strand -> Primer R 3' [1][3]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerf][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerf][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer F -> F Strand -> Primer RRC 3' [1][4]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerf][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerf][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerf][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer FRC -> F Strand -> FRC 3' [2][0]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerfrc][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerfrc][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer FRC -> F Strand -> Primer F 3' [2][1]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerfrc][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerfrc][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer FRC -> F Strand -> Primer FRC 3' [2][2]
                //5' and 3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerfrc][primerfrc], MPFR_RNDN);
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerfrc][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer FRC -> F Strand -> Primer R 3' [2][3]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerfrc][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerfrc][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer FRC -> F Strand -> Primer RRC 3' [2][4]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerfrc][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerfrc][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerfrc][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer R -> F Strand -> FRC 3' [3][0]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerr][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerr][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer R -> F Strand -> Primer F 3' [3][1]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerr][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerr][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer R -> F Strand -> Primer FRC 3' [3][2]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerr][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerr][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer R -> F Strand -> Primer R 3' [3][3]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerr][primerr], MPFR_RNDN);
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerr][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer R -> F Strand -> Primer RRC 3' [3][4]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerr][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerr][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerr][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer RRC -> F Strand -> FRC 3' [4][0]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerrrc][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerrrc][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][addr_end], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer RRC -> F Strand -> Primer F 3' [4][1]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerrrc][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerrrc][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][primerf], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer RRC -> F Strand -> Primer FRC 3' [4][2]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerrrc][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerrrc][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][primerfrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer RRC -> F Strand -> Primer R 3' [4][3]
                //5' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerrrc][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                //3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerrrc][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][primerr], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);

                //5' Primer RRC -> F Strand -> Primer RRC 3' [4][4]
                //5' and 3' Binding
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerrrc][primerrrc], MPFR_RNDN);
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].fstrand[primerrrc][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.tmp[3], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_mul_d(eq.tmp[3], eq.address_k_conc_vec[i].fstrand[primerrrc][primerrrc], eq_const, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
            }
            //Set k[K_FX] and k[K_RX] to average eq constant
            mpfr_div(eq.k[K_FX], eq.tmp[1], eq.tmp[0], MPFR_RNDN);
            mpfr_div(eq.k[K_RX], eq.tmp[2], eq.tmp[0], MPFR_RNDN);
            //All k values now set

            //Set initial concentrations
            //c0[F] and c0[R] are set by previous iteration or initialized before the loop
            //eq.tmp[0] still holds the total nonspecific concentration
            mpfr_set(eq.c0[X], eq.tmp[0], MPFR_RNDN);
            //All initial concs now set
            
            eq.solve_eq();

            //Now solve individual nonspecific concentrations and update c0 for next cycle
            for(unsigned int i = 0; i < addresses.size(); i++){
                mpfr_set_d(eq.address_k_conc_vec[i].last_f_conc, 0.0, MPFR_RNDN);
                mpfr_set_d(eq.address_k_conc_vec[i].last_r_conc, 0.0, MPFR_RNDN);

                for(int ii = 0; ii < 5; ii++){
                    for(int jj = 0; jj < 5; jj++){
                        mpfr_add(eq.address_k_conc_vec[i].last_f_conc, eq.address_k_conc_vec[i].last_f_conc, eq.address_k_conc_vec[i].fstrand[ii][jj], MPFR_RNDN);
                        mpfr_add(eq.address_k_conc_vec[i].last_r_conc, eq.address_k_conc_vec[i].last_r_conc, eq.address_k_conc_vec[i].rstrand[ii][jj], MPFR_RNDN);

                        mpfr_set_d(eq.address_k_conc_vec[i].fstrand_change[ii][jj], 0.0, MPFR_RNDN);
                        mpfr_set_d(eq.address_k_conc_vec[i].rstrand_change[ii][jj], 0.0, MPFR_RNDN);
                    }
                }
                //R Strand Bindings
                //[0][0]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[addr_end][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][addr_end], eq.address_k_conc_vec[i].fstrand_change[primerf][addr_end], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][addr_end], eq.address_k_conc_vec[i].fstrand_change[primerr][addr_end], eq.c[RX], MPFR_RNDN);
                
                //[0][1]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[addr_end][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][addr_end], eq.address_k_conc_vec[i].fstrand_change[primerf][addr_end], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][addr_end], eq.address_k_conc_vec[i].fstrand_change[primerr][addr_end], eq.c[RX], MPFR_RNDN);

                //[0][2]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[addr_end][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][addr_end], eq.address_k_conc_vec[i].fstrand_change[primerf][addr_end], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][addr_end], eq.address_k_conc_vec[i].fstrand_change[primerr][addr_end], eq.c[RX], MPFR_RNDN);

                //[0][3]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[addr_end][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][addr_end], eq.address_k_conc_vec[i].fstrand_change[primerf][addr_end], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][addr_end], eq.address_k_conc_vec[i].fstrand_change[primerr][addr_end], eq.c[RX], MPFR_RNDN);

                //[0][4]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[addr_end][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][addr_end], eq.address_k_conc_vec[i].fstrand_change[primerf][addr_end], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][addr_end], eq.address_k_conc_vec[i].fstrand_change[primerr][addr_end], eq.c[RX], MPFR_RNDN);

                //[1][0]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[primerf][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][primerfrc], eq.address_k_conc_vec[i].fstrand_change[primerf][primerfrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][primerfrc], eq.address_k_conc_vec[i].fstrand_change[primerr][primerfrc], eq.c[RX], MPFR_RNDN);

                //[1][1]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[primerf][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][primerfrc], eq.address_k_conc_vec[i].fstrand_change[primerf][primerfrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][primerfrc], eq.address_k_conc_vec[i].fstrand_change[primerr][primerfrc], eq.c[RX], MPFR_RNDN);

                //[1][2]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[primerf][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][primerfrc], eq.address_k_conc_vec[i].fstrand_change[primerf][primerfrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][primerfrc], eq.address_k_conc_vec[i].fstrand_change[primerr][primerfrc], eq.c[RX], MPFR_RNDN);

                //[1][3]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[primerf][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][primerfrc], eq.address_k_conc_vec[i].fstrand_change[primerf][primerfrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][primerfrc], eq.address_k_conc_vec[i].fstrand_change[primerr][primerfrc], eq.c[RX], MPFR_RNDN);

                //[1][4]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[primerf][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][primerfrc], eq.address_k_conc_vec[i].fstrand_change[primerf][primerfrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][primerfrc], eq.address_k_conc_vec[i].fstrand_change[primerr][primerfrc], eq.c[RX], MPFR_RNDN);

                //[2][0]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[primerfrc][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][primerf], eq.address_k_conc_vec[i].fstrand_change[primerf][primerf], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][primerf], eq.address_k_conc_vec[i].fstrand_change[primerr][primerf], eq.c[RX], MPFR_RNDN);

                //[2][1]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[primerfrc][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][primerf], eq.address_k_conc_vec[i].fstrand_change[primerf][primerf], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][primerf], eq.address_k_conc_vec[i].fstrand_change[primerr][primerf], eq.c[RX], MPFR_RNDN);

                //[2][3]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[primerfrc][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][primerf], eq.address_k_conc_vec[i].fstrand_change[primerf][primerf], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][primerf], eq.address_k_conc_vec[i].fstrand_change[primerr][primerf], eq.c[RX], MPFR_RNDN);

                //[3][0]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[primerr][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][primerrrc], eq.address_k_conc_vec[i].fstrand_change[primerf][primerrrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][primerrrc], eq.address_k_conc_vec[i].fstrand_change[primerr][primerrrc], eq.c[RX], MPFR_RNDN);

                //[3][1]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[primerr][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][primerrrc], eq.address_k_conc_vec[i].fstrand_change[primerf][primerrrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][primerrrc], eq.address_k_conc_vec[i].fstrand_change[primerr][primerrrc], eq.c[RX], MPFR_RNDN);

                //[3][2]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[primerr][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][primerrrc], eq.address_k_conc_vec[i].fstrand_change[primerf][primerrrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][primerrrc], eq.address_k_conc_vec[i].fstrand_change[primerr][primerrrc], eq.c[RX], MPFR_RNDN);

                //[3][3]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[primerr][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][primerrrc], eq.address_k_conc_vec[i].fstrand_change[primerf][primerrrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][primerrrc], eq.address_k_conc_vec[i].fstrand_change[primerr][primerrrc], eq.c[RX], MPFR_RNDN);

                //[3][4]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[primerr][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][primerrrc], eq.address_k_conc_vec[i].fstrand_change[primerf][primerrrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][primerrrc], eq.address_k_conc_vec[i].fstrand_change[primerr][primerrrc], eq.c[RX], MPFR_RNDN);

                //[4][0]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[primerrrc][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][primerr], eq.address_k_conc_vec[i].fstrand_change[primerf][primerr], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][primerr], eq.address_k_conc_vec[i].fstrand_change[primerr][primerr], eq.c[RX], MPFR_RNDN);

                //[4][1]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[primerrrc][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][primerr], eq.address_k_conc_vec[i].fstrand_change[primerf][primerr], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][primerr], eq.address_k_conc_vec[i].fstrand_change[primerr][primerr], eq.c[RX], MPFR_RNDN);

                //[4][3]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].rstrand[primerrrc][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerf][primerr], eq.address_k_conc_vec[i].fstrand_change[primerf][primerr], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].fstrand_change[primerr][primerr], eq.address_k_conc_vec[i].fstrand_change[primerr][primerr], eq.c[RX], MPFR_RNDN);


                //F Strand Bindings
                //[0][0]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[addr_end][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][addr_end], eq.address_k_conc_vec[i].rstrand_change[primerf][addr_end], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][addr_end], eq.address_k_conc_vec[i].rstrand_change[primerr][addr_end], eq.c[RX], MPFR_RNDN);
                
                //[0][1]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[addr_end][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][addr_end], eq.address_k_conc_vec[i].rstrand_change[primerf][addr_end], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][addr_end], eq.address_k_conc_vec[i].rstrand_change[primerr][addr_end], eq.c[RX], MPFR_RNDN);

                //[0][2]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[addr_end][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][addr_end], eq.address_k_conc_vec[i].rstrand_change[primerf][addr_end], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][addr_end], eq.address_k_conc_vec[i].rstrand_change[primerr][addr_end], eq.c[RX], MPFR_RNDN);

                //[0][3]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[addr_end][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][addr_end], eq.address_k_conc_vec[i].rstrand_change[primerf][addr_end], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][addr_end], eq.address_k_conc_vec[i].rstrand_change[primerr][addr_end], eq.c[RX], MPFR_RNDN);

                //[0][4]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[addr_end][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][addr_end], eq.address_k_conc_vec[i].rstrand_change[primerf][addr_end], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][addr_end], eq.address_k_conc_vec[i].rstrand_change[primerr][addr_end], eq.c[RX], MPFR_RNDN);

                //[1][0]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[primerf][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][primerfrc], eq.address_k_conc_vec[i].rstrand_change[primerf][primerfrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][primerfrc], eq.address_k_conc_vec[i].rstrand_change[primerr][primerfrc], eq.c[RX], MPFR_RNDN);

                //[1][1]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[primerf][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][primerfrc], eq.address_k_conc_vec[i].rstrand_change[primerf][primerfrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][primerfrc], eq.address_k_conc_vec[i].rstrand_change[primerr][primerfrc], eq.c[RX], MPFR_RNDN);

                //[1][2]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[primerf][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][primerfrc], eq.address_k_conc_vec[i].rstrand_change[primerf][primerfrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][primerfrc], eq.address_k_conc_vec[i].rstrand_change[primerr][primerfrc], eq.c[RX], MPFR_RNDN);

                //[1][3]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[primerf][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][primerfrc], eq.address_k_conc_vec[i].rstrand_change[primerf][primerfrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][primerfrc], eq.address_k_conc_vec[i].rstrand_change[primerr][primerfrc], eq.c[RX], MPFR_RNDN);

                //[1][4]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[primerf][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][primerfrc], eq.address_k_conc_vec[i].rstrand_change[primerf][primerfrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][primerfrc], eq.address_k_conc_vec[i].rstrand_change[primerr][primerfrc], eq.c[RX], MPFR_RNDN);

                //[2][0]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[primerfrc][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][primerf], eq.address_k_conc_vec[i].rstrand_change[primerf][primerf], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][primerf], eq.address_k_conc_vec[i].rstrand_change[primerr][primerf], eq.c[RX], MPFR_RNDN);

                //[2][1]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[primerfrc][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][primerf], eq.address_k_conc_vec[i].rstrand_change[primerf][primerf], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][primerf], eq.address_k_conc_vec[i].rstrand_change[primerr][primerf], eq.c[RX], MPFR_RNDN);

                //[2][3]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[primerfrc][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][primerf], eq.address_k_conc_vec[i].rstrand_change[primerf][primerf], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][primerf], eq.address_k_conc_vec[i].rstrand_change[primerr][primerf], eq.c[RX], MPFR_RNDN);

                //[3][0]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[primerr][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][primerrrc], eq.address_k_conc_vec[i].rstrand_change[primerf][primerrrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][primerrrc], eq.address_k_conc_vec[i].rstrand_change[primerr][primerrrc], eq.c[RX], MPFR_RNDN);

                //[3][1]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[primerr][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][primerrrc], eq.address_k_conc_vec[i].rstrand_change[primerf][primerrrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][primerrrc], eq.address_k_conc_vec[i].rstrand_change[primerr][primerrrc], eq.c[RX], MPFR_RNDN);

                //[3][2]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[primerr][primerfrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_frc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][primerrrc], eq.address_k_conc_vec[i].rstrand_change[primerf][primerrrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][primerrrc], eq.address_k_conc_vec[i].rstrand_change[primerr][primerrrc], eq.c[RX], MPFR_RNDN);

                //[3][3]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[primerr][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][primerrrc], eq.address_k_conc_vec[i].rstrand_change[primerf][primerrrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][primerrrc], eq.address_k_conc_vec[i].rstrand_change[primerr][primerrrc], eq.c[RX], MPFR_RNDN);

                //[3][4]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[primerr][primerrrc], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][primerrrc], eq.address_k_conc_vec[i].rstrand_change[primerf][primerrrc], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][primerrrc], eq.address_k_conc_vec[i].rstrand_change[primerr][primerrrc], eq.c[RX], MPFR_RNDN);

                //[4][0]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[primerrrc][addr_end], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_f_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[i].dhds_primer_r_addr_rrc, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][primerr], eq.address_k_conc_vec[i].rstrand_change[primerf][primerr], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][primerr], eq.address_k_conc_vec[i].rstrand_change[primerr][primerr], eq.c[RX], MPFR_RNDN);

                //[4][1]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[primerrrc][primerf], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_f, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][primerr], eq.address_k_conc_vec[i].rstrand_change[primerf][primerr], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][primerr], eq.address_k_conc_vec[i].rstrand_change[primerr][primerr], eq.c[RX], MPFR_RNDN);

                //[4][3]
                mpfr_set(eq.c0[X], eq.address_k_conc_vec[i].fstrand[primerrrc][primerr], MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_f_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_FX], eq_const, MPFR_RNDN);
                eq_const = dhds_to_eq_const(eq.address_k_conc_vec[addr].dhds_primer_r_addr_r, temp_c_profile[cycle-1]);
                mpfr_set_d(eq.k[K_RX], eq_const, MPFR_RNDN);
                eq.calc_cx();
                eq.calc_bound_concs();
                //Remove bound primer from primer concnetration
                mpfr_sub(eq.c0[F], eq.c0[F], eq.c[FX], MPFR_RNDN);
                mpfr_sub(eq.c0[R], eq.c0[R], eq.c[RX], MPFR_RNDN);
                //Add to concentration of elongated strand
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerf][primerr], eq.address_k_conc_vec[i].rstrand_change[primerf][primerr], eq.c[FX], MPFR_RNDN);
                mpfr_add(eq.address_k_conc_vec[i].rstrand_change[primerr][primerr], eq.address_k_conc_vec[i].rstrand_change[primerr][primerr], eq.c[RX], MPFR_RNDN);

                //Add changes to strand concentrations and sum all f and r strands
                mpfr_set_d(eq.address_k_conc_vec[i].total_f_conc, 0.0, MPFR_RNDN);
                mpfr_set_d(eq.address_k_conc_vec[i].total_r_conc, 0.0, MPFR_RNDN);
                for(int ii = 0; ii < 5; ii++){
                    for(int jj = 0; jj < 5; jj++){
                        mpfr_add(eq.address_k_conc_vec[i].fstrand[ii][jj], eq.address_k_conc_vec[i].fstrand[ii][jj], eq.address_k_conc_vec[i].fstrand_change[ii][jj], MPFR_RNDN);
                        mpfr_add(eq.address_k_conc_vec[i].rstrand[ii][jj], eq.address_k_conc_vec[i].rstrand[ii][jj], eq.address_k_conc_vec[i].rstrand_change[ii][jj], MPFR_RNDN);
                        mpfr_add(eq.address_k_conc_vec[i].total_f_conc, eq.address_k_conc_vec[i].total_f_conc, eq.address_k_conc_vec[i].fstrand[ii][jj], MPFR_RNDN);
                        mpfr_add(eq.address_k_conc_vec[i].total_r_conc, eq.address_k_conc_vec[i].total_r_conc, eq.address_k_conc_vec[i].rstrand[ii][jj], MPFR_RNDN);
                    }
                }

            }
            

            /*
            //tmp[0] = percent bound forward primer binding sites
            mpfr_add(eq.tmp[3], eq.c[FA], eq.c[RA], MPFR_RNDN);
            mpfr_set(eq.tmp[2], eq.tmp[3], MPFR_RNDN);
            mpfr_add(eq.tmp[3], eq.tmp[3], eq.c[A], MPFR_RNDN);
            mpfr_div(eq.tmp[0], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
            mpfr_mul_d(eq.tmp[0], eq.tmp[0], 100., MPFR_RNDN);

            //tmp[1] = percent bound rev primer binding sites
            mpfr_add(eq.tmp[3], eq.c[FB], eq.c[RB], MPFR_RNDN);
            mpfr_set(eq.tmp[2], eq.tmp[3], MPFR_RNDN);
            mpfr_add(eq.tmp[3], eq.tmp[3], eq.c[B], MPFR_RNDN);
            mpfr_div(eq.tmp[1], eq.tmp[2], eq.tmp[3], MPFR_RNDN);
            mpfr_mul_d(eq.tmp[1], eq.tmp[1], 100., MPFR_RNDN);
            */
            mpfr_div(eq.tmp[0], eq.address_k_conc_vec[addr].total_f_conc, eq.address_k_conc_vec[addr].last_f_conc, MPFR_RNDN); 
            mpfr_sub_d(eq.tmp[0], eq.tmp[0], 1.0, MPFR_RNDN);

            mpfr_div(eq.tmp[1], eq.address_k_conc_vec[addr].total_r_conc, eq.address_k_conc_vec[addr].last_r_conc, MPFR_RNDN); 
            mpfr_sub_d(eq.tmp[1], eq.tmp[1], 1.0, MPFR_RNDN);

            FILE *outfile = fopen(out_filename, "a");
            //cycle,f_primer,r_primer,spec_f,spec_r,nonspec_f,nonspec_r
            mpfr_fprintf(outfile, "%02u,%lf,%.9Rf,%.9Rf,%.9Re,%.9Re,%.9Re,%.9Re", cycle, temp_c_profile[cycle-1], eq.tmp[0], eq.tmp[1], eq.c0[F], eq.c0[R], eq.address_k_conc_vec[addr].total_f_conc, eq.address_k_conc_vec[addr].total_r_conc);
            mpfr_set_d(eq.tmp[0], 0., MPFR_RNDN);
            mpfr_set_d(eq.tmp[1], 0., MPFR_RNDN);
            mpfr_set_d(eq.tmp[2], 0., MPFR_RNDN);
            mpfr_set_d(eq.tmp[3], 0., MPFR_RNDN);
            for(unsigned int i = 0; i < addresses.size(); i++){
                if (i == addr)
                    continue;
                mpfr_add(eq.tmp[0], eq.tmp[0], eq.address_k_conc_vec[i].total_f_conc, MPFR_RNDN);
                mpfr_add(eq.tmp[1], eq.tmp[1], eq.address_k_conc_vec[i].last_f_conc, MPFR_RNDN);
                mpfr_add(eq.tmp[2], eq.tmp[2], eq.address_k_conc_vec[i].total_r_conc, MPFR_RNDN);
                mpfr_add(eq.tmp[3], eq.tmp[3], eq.address_k_conc_vec[i].last_r_conc, MPFR_RNDN);
            }
            /*
            mpfr_add_d(eq.avg_nonspec_amp, eq.avg_nonspec_amp, 1.0, MPFR_RNDN);
            mpfr_mul(eq.tmp[2], eq.tmp[0], eq.avg_nonspec_amp, MPFR_RNDN);
            mpfr_mul(eq.tmp[3], eq.tmp[1], eq.avg_nonspec_amp, MPFR_RNDN);
            */
            mpfr_div(eq.tmp[1], eq.tmp[0], eq.tmp[1],MPFR_RNDN);
            mpfr_sub_d(eq.tmp[1], eq.tmp[1], 1., MPFR_RNDN);
            mpfr_div(eq.tmp[3], eq.tmp[2], eq.tmp[3],MPFR_RNDN);
            mpfr_sub_d(eq.tmp[3], eq.tmp[3], 1., MPFR_RNDN);
            mpfr_set(eq.last_nonspec_frc_total, eq.tmp[0], MPFR_RNDN);
            mpfr_set(eq.last_nonspec_rrc_total, eq.tmp[1], MPFR_RNDN);
            mpfr_fprintf(outfile, ",%.9Re,%.9Rf,%.9Re,%.9Rf\n", eq.tmp[0], eq.tmp[1], eq.tmp[2], eq.tmp[3]);
            fclose(outfile);
        }
    }

}
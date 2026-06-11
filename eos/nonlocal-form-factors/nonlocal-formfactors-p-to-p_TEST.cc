/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019 Nico Gubernari
 * Copyright (c) 2021 Méril Reboud
 * Copyright (c) 2025 Danny van Dyk
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <test/test.hh>
#include <eos/nonlocal-form-factors/nonlocal-formfactors.hh>

#include <iomanip>
#include <limits>

using namespace test;
using namespace eos;

class NonlocalFormFactorGvDV2020Test :
    public TestCase
{
    public:
        NonlocalFormFactorGvDV2020Test() :
            TestCase("nonlocal_formfactor_GvDV2020_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["mass::B_d"]                               = 5.279;
                p["mass::K_d"]                               = 0.492;
                p["mass::J/psi"]                             = 3.0969;
                p["mass::psi(2S)"]                           = 3.6860;
                p["mass::D^0"]                               = 1.86723;
                p["b->sccbar::t_0"]                          = 4.0;
                p["b->sccbar::t_s"]                          = -17.4724;
                p["b->sccbar::chiOPE@GvDV2020"]              = 1.81e-4;
                p["B->Kccbar::Re{alpha_0^plus}@GvDV2020"]    = 2.0;
                p["B->Kccbar::Im{alpha_0^plus}@GvDV2020"]    = 3.0;
                p["B->Kccbar::Re{alpha_1^plus}@GvDV2020"]    = 4.0;
                p["B->Kccbar::Im{alpha_1^plus}@GvDV2020"]    = 5.0;
                p["B->Kccbar::Re{alpha_2^plus}@GvDV2020"]    = 6.0;
                p["B->Kccbar::Im{alpha_2^plus}@GvDV2020"]    = 7.0;

                Options o = { { "model"_ok, "WET"_ov } };

                auto nff = NonlocalFormFactor<PToP>::make("B->K::GvDV2020", p, o);


                auto diagnostics = nff->diagnostics();

                std::cout << "Diagnostics:" << std::endl;
                for (auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                std::cout << "Diagnostics ended" << std::endl;

                static const std::vector<std::pair<double, double>> reference
                {
                    /* outer functions */
                    std::make_pair(   0.0,      eps),            // Re{1/phi_+(q2 = 0.0)}
                    std::make_pair(   0.0,      eps),            // Im{1/phi_+(q2 = 0.0)}
                    std::make_pair( -17.44509,  eps),            // Re{phi_+(q2 = 16.0)}
                    std::make_pair(   4.863096, eps),            // Im{phi_+(q2 = 16.0)}

                    std::make_pair( -18.00857,  eps),            // Re{PGvDV2020(q2 = 1.0, sXY = 0.6+0.8i, {2.0, 3.0, 4.0, 5.0})}
                    std::make_pair(   0.0,      eps),            // Im{PGvDV2020(q2 = 1.0, sXY = 0.6+0.8i, {2.0, 3.0, 4.0, 5.0})}


                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_RELATIVE_ERROR(real(nff->H_plus(-1.0)),  0.0925106365,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_plus(-1.0)),  0.1116021529,   eps);
                TEST_CHECK_NEARLY_EQUAL(  real(nff->H_plus(0.0)),   0.,             eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(nff->H_plus(0.0)),   0.,             eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_plus(4.0)),  -0.7306769592,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_plus(4.0)),  -0.8733733008,   eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_plus(12.0)),  7.9949120618,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_plus(12.0)),  9.3613667608,   eps);

                TEST_CHECK_RELATIVE_ERROR(real(nff->H_plus_residue_jpsi()),   11.5294356,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_plus_residue_jpsi()),   13.5982909,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_plus_residue_psi2s()), -3.10770287,  eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_plus_residue_psi2s()), -3.61681383,  eps);
            }

            {
                Parameters p = Parameters::Defaults();
                p["mass::B_d"]                                   =  5.279;
                p["mass::K_d"]                                   =  0.492;
                p["mass::J/psi"]                                 =  3.0969;
                p["mass::psi(2S)"]                               =  3.6860;
                p["mass::D^0"]                                   =  1.86723;
                p["b->sccbar::t_0"]                              =  4.0;
                p["b->sccbar::t_s"]                              = -17.4724;
                p["b->sccbar::chiOPE@GvDV2020"]                  =  1.81e-4;
                p["B->Kccbar::Re_Hhat_at_m7_plus@GRvDV2022"]     =  0.01;
                p["B->Kccbar::Im_Hhat_at_m7_plus@GRvDV2022"]     = -0.02;
                p["B->Kccbar::Re_Hhat_at_m5_plus@GRvDV2022"]     = -0.03;
                p["B->Kccbar::Im_Hhat_at_m5_plus@GRvDV2022"]     =  0.04;
                p["B->Kccbar::Re_Hhat_at_m3_plus@GRvDV2022"]     = -0.05;
                p["B->Kccbar::Im_Hhat_at_m3_plus@GRvDV2022"]     = -0.06;
                p["B->Kccbar::Re_Hhat_at_m1_plus@GRvDV2022"]     =  0.07;
                p["B->Kccbar::Im_Hhat_at_m1_plus@GRvDV2022"]     = -0.08;
                p["B->Kccbar::Abs_Hhat_at_Jpsi_plus@GRvDV2022"]  =  0.09;
                p["B->Kccbar::Arg_Hhat_at_Jpsi_plus@GRvDV2022"]  =  0.10;
                p["B->Kccbar::Abs_Hhat_at_psi2S_plus@GRvDV2022"] =  0.11;
                p["B->Kccbar::Arg_Hhat_at_psi2S_plus@GRvDV2022"] = -0.12;


                Options o = { { "model"_ok, "WET"_ov } };

                auto nff = NonlocalFormFactor<PToP>::make("B->K::GRvDV2022order5", p, o);

                TEST_CHECK_RELATIVE_ERROR(real(nff->H_plus(-1.0)),  0.00124139, eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_plus(-1.0)), -0.00141873, eps);
                TEST_CHECK_NEARLY_EQUAL(  real(nff->H_plus(0.0)),   0.,         eps);
                TEST_CHECK_NEARLY_EQUAL(  imag(nff->H_plus(0.0)),   0.,         eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_plus(4.0)),  -0.221699,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_plus(4.0)),  -0.264192,   eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_plus(12.0)), -20.9153,    eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_plus(12.0)), -37.2671,    eps);

                TEST_CHECK_RELATIVE_ERROR(real(nff->H_plus_residue_jpsi()),   0.122788,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_plus_residue_jpsi()),   0.0123199,  eps);
                TEST_CHECK_RELATIVE_ERROR(real(nff->H_plus_residue_psi2s()), -0.019915,   eps);
                TEST_CHECK_RELATIVE_ERROR(imag(nff->H_plus_residue_psi2s()),  0.00240133, eps);

                TEST_CHECK_RELATIVE_ERROR(nff->weak_bound(),   149077281.62, eps);
                TEST_CHECK_RELATIVE_ERROR(nff->strong_bound(), 240635402.59, eps);
            }

        }
} nonlocal_formfactor_gvdv2020_test;



class NonlocalFormFactorGRvDV2022order5Test :
    public TestCase
{
    public:
        NonlocalFormFactorGRvDV2022order5Test() :
            TestCase("nonlocal_formfactor_GRvDV2022order5_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["mass::B_d"]                               = 5.279;
                p["mass::K_d"]                               = 0.492;
                p["mass::J/psi"]                             = 3.0969;
                p["mass::psi(2S)"]                           = 3.6860;
                p["mass::D^0"]                               = 1.86723;
                p["b->sccbar::t_0"]                          = 4.0;
                p["b->sccbar::t_s"]                          = -17.4724;
                p["b->sccbar::chiOPE@GvDV2020"]        = 1.81e-4;
                p["B->Kccbar::Re_Hhat_at_m7_plus@GRvDV2022"]     =  0.01;
                p["B->Kccbar::Im_Hhat_at_m7_plus@GRvDV2022"]     = -0.02;
                p["B->Kccbar::Re_Hhat_at_m5_plus@GRvDV2022"]     = -0.03;
                p["B->Kccbar::Im_Hhat_at_m5_plus@GRvDV2022"]     =  0.04;
                p["B->Kccbar::Re_Hhat_at_m3_plus@GRvDV2022"]     = -0.05;
                p["B->Kccbar::Im_Hhat_at_m3_plus@GRvDV2022"]     = -0.06;
                p["B->Kccbar::Re_Hhat_at_m1_plus@GRvDV2022"]     =  0.07;
                p["B->Kccbar::Im_Hhat_at_m1_plus@GRvDV2022"]     = -0.08;
                p["B->Kccbar::Abs_Hhat_at_Jpsi_plus@GRvDV2022"]  =  0.09;
                p["B->Kccbar::Arg_Hhat_at_Jpsi_plus@GRvDV2022"]  =  0.10;
                p["B->Kccbar::Abs_Hhat_at_psi2S_plus@GRvDV2022"] =  0.11;
                p["B->Kccbar::Arg_Hhat_at_psi2S_plus@GRvDV2022"] = -0.12;

                Options o = { { "model"_ok, "WET" } };

                auto nff = NonlocalFormFactor<PToP>::make("B->K::GRvDV2022order5", p, o);


                auto diagnostics = nff->diagnostics();

                std::cout << "Diagnostics:" << std::endl;
                for (auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                std::cout << "Diagnostics ended" << std::endl;

                static const std::vector<std::pair<double, double>> reference
                {
                    /* outer functions */
                    std::make_pair(  5.14673, eps),             // Re{phi_+(q2 = 0.0)}
                    std::make_pair(  0.0, eps),             // Im{phi_+(q2 = 0.0)}
                    std::make_pair(  5.43688,  eps),             // Re{phi_+(q2 = -4.0)}
                    std::make_pair(  0.0,  eps),             // Im{phi_+(q2 = -4.0)}
                    std::make_pair(  6.58341,  eps),             // Re{phi_+(q2 = 7.0)}
                    std::make_pair(  0.0,  eps),             // Im{phi_+(q2 = 7.0)}
                    std::make_pair(  1.67205,  eps),             // Re{phi_+(q2 = 16.0)}
                    std::make_pair(  25.5052,  eps),             // Im{phi_+(q2 = 16.0)}

                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                // TEST_CHECK_RELATIVE_ERROR(real(nff->H_plus(-1.0)),  0.0925106365,   eps);
                // TEST_CHECK_RELATIVE_ERROR(imag(nff->H_plus(-1.0)),  0.1116021529,   eps);
                // TEST_CHECK_NEARLY_EQUAL(  real(nff->H_plus(0.0)),   0.,             eps);
                // TEST_CHECK_NEARLY_EQUAL(  imag(nff->H_plus(0.0)),   0.,             eps);
                // TEST_CHECK_RELATIVE_ERROR(real(nff->H_plus(4.0)),  -0.7306769592,   eps);
                // TEST_CHECK_RELATIVE_ERROR(imag(nff->H_plus(4.0)),  -0.8733733008,   eps);
                // TEST_CHECK_RELATIVE_ERROR(real(nff->H_plus(12.0)),  7.9949120618,   eps);
                // TEST_CHECK_RELATIVE_ERROR(imag(nff->H_plus(12.0)),  9.3613667608,   eps);

                // TEST_CHECK_RELATIVE_ERROR(real(nff->H_plus_residue_jpsi()),   11.5294356,  eps);
                // TEST_CHECK_RELATIVE_ERROR(imag(nff->H_plus_residue_jpsi()),   13.5982909,  eps);
                // TEST_CHECK_RELATIVE_ERROR(real(nff->H_plus_residue_psi2s()), -3.10770287,  eps);
                // TEST_CHECK_RELATIVE_ERROR(imag(nff->H_plus_residue_psi2s()), -3.61681383,  eps);
            }
        }
} nonlocal_formfactor_GRvDV2022order5_test;



class NonlocalFormFactorGRV2026Test :
    public TestCase
{
    public:
        NonlocalFormFactorGRV2026Test() :
            TestCase("nonlocal_formfactor_GRV2026_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["mass::B_d"]                               = 5.279;
                p["mass::K_d"]                               = 0.492;
                p["mass::J/psi"]                             = 3.0969;
                p["mass::psi(2S)"]                           = 3.6860;
                p["mass::D^0"]                               = 1.86723;
                p["b->sccbar::t_0"]                          = 4.0;
                p["b->sccbar::t_s"]                          = -17.4724;
                p["b->sccbar::t_V@GRV2026"]                  = 13.910363715599999; // DDbar subthreshold
                p["b->sccbar::chi_tilde_OPE_V@GRV2026"]        = +5.796e-5;
                p["B->Kccbar::Re_Hhat_at_m7_plus@GRvDV2022"]     =  0.01;
                p["B->Kccbar::Im_Hhat_at_m7_plus@GRvDV2022"]     = -0.02;
                p["B->Kccbar::Re_Hhat_at_m5_plus@GRvDV2022"]     = -0.03;
                p["B->Kccbar::Im_Hhat_at_m5_plus@GRvDV2022"]     =  0.04;
                p["B->Kccbar::Re_Hhat_at_m3_plus@GRvDV2022"]     = -0.05;
                p["B->Kccbar::Im_Hhat_at_m3_plus@GRvDV2022"]     = -0.06;
                p["B->Kccbar::Re_Hhat_at_m1_plus@GRvDV2022"]     =  0.07;
                p["B->Kccbar::Im_Hhat_at_m1_plus@GRvDV2022"]     = -0.08;
                p["B->Kccbar::Abs_Hhat_at_Jpsi_plus@GRvDV2022"]  =  0.09;
                p["B->Kccbar::Arg_Hhat_at_Jpsi_plus@GRvDV2022"]  =  0.10;
                p["B->Kccbar::Re_Hhat_at_t0_plus@GRvDV2022"] =  0.11;
                p["B->Kccbar::Im_Hhat_at_t0_plus@GRvDV2022"] = -0.12;

                Options o = { { "model"_ok, "WET" } };

                auto nff = NonlocalFormFactor<PToP>::make("B->K::GRV2026", p, o);


                auto diagnostics = nff->diagnostics();

                std::cout << std::setprecision(17);
                std::cout << "Diagnostics:" << std::endl;
                for (auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                std::cout << "Diagnostics ended" << std::endl;

                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair( -17.4724 , eps ),		// t_s
					std::make_pair( 13.910363715599999 , eps ),		// s_G
					std::make_pair( -0.20467397269356458 , eps ),		// Re_z_poles
					std::make_pair( 0 , eps ),		// Im_z_poles
					std::make_pair( -0.69384244249996163 , eps ),		// Re_z_poles
					std::make_pair( 0 , eps ),		// Im_z_poles
					std::make_pair( 5.4154 , eps ),		// Mres
					std::make_pair( -61.874744727032436 , eps ),		// Re{phi_+(q2 = -7.0)}
					std::make_pair( -7.5774708078019861e-15 , eps ),		// Im{phi_+(q2 = -7.0)}
					std::make_pair( -388.41456644867367 , eps ),		// Re{phi_+(q2 = -1.0)}
					std::make_pair( -4.7567065554357515e-14 , eps ),		// Im{phi_+(q2 = -1.0)}
					std::make_pair( 69.429346919727379 , eps ),		// Re{phi_+(q2 = 5.0)}
					std::make_pair( 0 , eps ),		// Im{phi_+(q2 = 5.0)}
					std::make_pair( 22.681788564918079 , eps ),		// Re{phi_+(q2 = 16.0)}
					std::make_pair( 22.723306968408291 , eps ),		// Im{phi_+(q2 = 16.0)}
					std::make_pair( 1.089969122226357 , eps ),		// Re{phi_+(q2 = 25.0)}
					std::make_pair( 0.2031583642442763 , eps ),		// Im{phi_+(q2 = 25.0)}
					std::make_pair( -0.21859161591241955 , eps ),		// Re{phi_+(q2 = 35.0)}
					std::make_pair( 0.67332270237146974 , eps ),		// Im{phi_+(q2 = 35.0)}
					std::make_pair( 0.01 , eps ),		// Re{P_GRV2026_lagrange(q2 = -7.000000)}
					std::make_pair( -0.02 , eps ),		// Im{P_GRV2026_lagrange(q2 = -7.000000)}
					std::make_pair( 0.010000000000001563 , eps ),		// Re{P_GRV2026_monomial(q2 = -7.000000)}
					std::make_pair( -0.019999999999999574 , eps ),		// Im{P_GRV2026_monomial(q2 = -7.000000)}
					std::make_pair( 0.070000000000000007 , eps ),		// Re{P_GRV2026_lagrange(q2 = -1.000000)}
					std::make_pair( -0.080000000000000002 , eps ),		// Im{P_GRV2026_lagrange(q2 = -1.000000)}
					std::make_pair( 0.070000000000000562 , eps ),		// Re{P_GRV2026_monomial(q2 = -1.000000)}
					std::make_pair( -0.080000000000000071 , eps ),		// Im{P_GRV2026_monomial(q2 = -1.000000)}
					std::make_pair( -0.56645857060104632 , eps ),		// Re{P_GRV2026_lagrange(q2 = 5.000000)}
					std::make_pair( -0.84030373192302132 , eps ),		// Im{P_GRV2026_lagrange(q2 = 5.000000)}
					std::make_pair( -0.56645857060104632 , eps ),		// Re{P_GRV2026_monomial(q2 = 5.000000)}
					std::make_pair( -0.84030373192302177 , eps ),		// Im{P_GRV2026_monomial(q2 = 5.000000)}
					std::make_pair( -50900.598014783165 , eps ),		// Re{P_GRV2026_lagrange(q2 = 16.000000)}
					std::make_pair( 4901.7099390189987 , eps ),		// Im{P_GRV2026_lagrange(q2 = 16.000000)}
					std::make_pair( -50900.598014783172 , eps ),		// Re{P_GRV2026_monomial(q2 = 16.000000)}
					std::make_pair( 4901.7099390190033 , eps ),		// Im{P_GRV2026_monomial(q2 = 16.000000)}
					std::make_pair( 36303.702108676618 , eps ),		// Re{P_GRV2026_lagrange(q2 = 25.000000)}
					std::make_pair( -27914.180154729282 , eps ),		// Im{P_GRV2026_lagrange(q2 = 25.000000)}
					std::make_pair( 36303.702108676604 , eps ),		// Re{P_GRV2026_monomial(q2 = 25.000000)}
					std::make_pair( -27914.180154729285 , eps ),		// Im{P_GRV2026_monomial(q2 = 25.000000)}
					std::make_pair( -24664.679150138709 , eps ),		// Re{P_GRV2026_lagrange(q2 = 35.000000)}
					std::make_pair( -34188.00268331109 , eps ),		// Im{P_GRV2026_lagrange(q2 = 35.000000)}
					std::make_pair( -24664.679150138705 , eps ),		// Re{P_GRV2026_monomial(q2 = 35.000000)}
					std::make_pair( -34188.002683311082 , eps ),		// Im{P_GRV2026_monomial(q2 = 35.000000)}
					std::make_pair( 0.11 , eps ),		// Re{GRV2026_monomial_coeff_0}
					std::make_pair( -0.12 , eps ),		// Im{GRV2026_monomial_coeff_0}
					std::make_pair( 17.678674364094149 , eps ),		// Re{GRV2026_monomial_coeff_1}
					std::make_pair( 18.083644293453311 , eps ),		// Im{GRV2026_monomial_coeff_1}
					std::make_pair( -281.42636186568086 , eps ),		// Re{GRV2026_monomial_coeff_2}
					std::make_pair( -317.97841867708956 , eps ),		// Im{GRV2026_monomial_coeff_2}
					std::make_pair( 596.86287287946664 , eps ),		// Re{GRV2026_monomial_coeff_3}
					std::make_pair( 1013.1900932274359 , eps ),		// Im{GRV2026_monomial_coeff_3}
					std::make_pair( 6660.204708579321 , eps ),		// Re{GRV2026_monomial_coeff_4}
					std::make_pair( 7535.9170900574536 , eps ),		// Im{GRV2026_monomial_coeff_4}
					std::make_pair( -24547.194862145938 , eps ),		// Re{GRV2026_monomial_coeff_5}
					std::make_pair( -35116.75852858369 , eps ),		// Im{GRV2026_monomial_coeff_5}
					std::make_pair( -0.00055342207322183856 , eps ),		// Re{H_+(q2 = -7.000000)}
					std::make_pair( 0.0011068441464436773 , eps ),		// Im{H_+(q2 = -7.000000)}
					std::make_pair( -0.00080787244873840717 , eps ),		// Re{H_+(q2 = -1.000000)}
					std::make_pair( 0.00092328279855817996 , eps ),		// Im{H_+(q2 = -1.000000)}
					std::make_pair( -0.067025484719303363 , eps ),		// Re{H_+(q2 = 5.000000)}
					std::make_pair( -0.099427862630482061 , eps ),		// Im{H_+(q2 = 5.000000)}
					std::make_pair( -644.89657639436928 , eps ),		// Re{H_+(q2 = 16.000000)}
					std::make_pair( -1456.3116495828838 , eps ),		// Im{H_+(q2 = 16.000000)}
					std::make_pair( -30960.574402632617 , eps ),		// Re{H_+(q2 = 25.000000)}
					std::make_pair( 27338.867982164866 , eps ),		// Im{H_+(q2 = 25.000000)}
					std::make_pair( 53237.34728653184 , eps ),		// Re{H_+(q2 = 35.000000)}
					std::make_pair( -26683.150686309429 , eps ),		// Im{H_+(q2 = 35.000000)}
					std::make_pair( 1835751505.1550193 , eps ),		// GRV2026 weak bound
					std::make_pair( 1938463628.4095731 , eps ),		// GRV2026 strong bound
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                // -------------------Different subthreshold-------------------

                p["b->sccbar::t_V@GRV2026"] = 10.445027650378242; // J/psi pi subthreshold

                auto nff2 = NonlocalFormFactor<PToP>::make("B->K::GRV2026", p, o);

                auto diagnostics2 = nff2->diagnostics();

                std::cout << "Diagnostics:" << std::endl;
                for (auto & d : diagnostics2)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                std::cout << "Diagnostics ended" << std::endl;

                static const std::vector<std::pair<double, double>> reference2
                {
                    /* outer functions */
                    std::make_pair( -17.4724 , eps ),		// t_s
					std::make_pair( 10.445027650378242 , eps ),		// s_G
					std::make_pair( -0.46620747262012496 , eps ),		// Re_z_poles
					std::make_pair( 0 , eps ),		// Im_z_poles
					std::make_pair( 3.6859999999999999 , eps ),		// Mres
					std::make_pair( 5.4154 , eps ),		// Mres
					std::make_pair( -23.809485784560046 , eps ),		// Re{phi_+(q2 = -7.0)}
					std::make_pair( -2.9158210555405866e-15 , eps ),		// Im{phi_+(q2 = -7.0)}
					std::make_pair( -135.77424625055733 , eps ),		// Re{phi_+(q2 = -1.0)}
					std::make_pair( -1.6627549607738956e-14 , eps ),		// Im{phi_+(q2 = -1.0)}
					std::make_pair( 21.300493463109653 , eps ),		// Re{phi_+(q2 = 5.0)}
					std::make_pair( 0 , eps ),		// Im{phi_+(q2 = 5.0)}
					std::make_pair( 2.9226271941942188 , eps ),		// Re{phi_+(q2 = 16.0)}
					std::make_pair( -1.5141421546023177 , eps ),		// Im{phi_+(q2 = 16.0)}
					std::make_pair( 0.14399945914555839 , eps ),		// Re{phi_+(q2 = 25.0)}
					std::make_pair( -0.32461685293221776 , eps ),		// Im{phi_+(q2 = 25.0)}
					std::make_pair( 0.26317861196563697 , eps ),		// Re{phi_+(q2 = 35.0)}
					std::make_pair( 0.20593257143494992 , eps ),		// Im{phi_+(q2 = 35.0)}
					std::make_pair( 0.01 , eps ),		// Re{P_GRV2026_lagrange(q2 = -7.000000)}
					std::make_pair( -0.02 , eps ),		// Im{P_GRV2026_lagrange(q2 = -7.000000)}
					std::make_pair( 0.010000000000004672 , eps ),		// Re{P_GRV2026_monomial(q2 = -7.000000)}
					std::make_pair( -0.02000000000000135 , eps ),		// Im{P_GRV2026_monomial(q2 = -7.000000)}
					std::make_pair( 0.070000000000000007 , eps ),		// Re{P_GRV2026_lagrange(q2 = -1.000000)}
					std::make_pair( -0.080000000000000002 , eps ),		// Im{P_GRV2026_lagrange(q2 = -1.000000)}
					std::make_pair( 0.070000000000000062 , eps ),		// Re{P_GRV2026_monomial(q2 = -1.000000)}
					std::make_pair( -0.080000000000000571 , eps ),		// Im{P_GRV2026_monomial(q2 = -1.000000)}
					std::make_pair( -1.1886868110991453 , eps ),		// Re{P_GRV2026_lagrange(q2 = 5.000000)}
					std::make_pair( -1.7766421049372632 , eps ),		// Im{P_GRV2026_lagrange(q2 = 5.000000)}
					std::make_pair( -1.1886868110991449 , eps ),		// Re{P_GRV2026_monomial(q2 = 5.000000)}
					std::make_pair( -1.7766421049372623 , eps ),		// Im{P_GRV2026_monomial(q2 = 5.000000)}
					std::make_pair( 11073.737514228895 , eps ),		// Re{P_GRV2026_lagrange(q2 = 16.000000)}
					std::make_pair( -531.32336494402841 , eps ),		// Im{P_GRV2026_lagrange(q2 = 16.000000)}
					std::make_pair( 11073.737514228893 , eps ),		// Re{P_GRV2026_monomial(q2 = 16.000000)}
					std::make_pair( -531.32336494402898 , eps ),		// Im{P_GRV2026_monomial(q2 = 16.000000)}
					std::make_pair( -6408.5108147018072 , eps ),		// Re{P_GRV2026_lagrange(q2 = 25.000000)}
					std::make_pair( -7524.9377030634842 , eps ),		// Im{P_GRV2026_lagrange(q2 = 25.000000)}
					std::make_pair( -6408.5108147018054 , eps ),		// Re{P_GRV2026_monomial(q2 = 25.000000)}
					std::make_pair( -7524.937703063486 , eps ),		// Im{P_GRV2026_monomial(q2 = 25.000000)}
					std::make_pair( -8812.8552414374608 , eps ),		// Re{P_GRV2026_lagrange(q2 = 35.000000)}
					std::make_pair( 2226.4290588318459 , eps ),		// Im{P_GRV2026_lagrange(q2 = 35.000000)}
					std::make_pair( -8812.8552414374608 , eps ),		// Re{P_GRV2026_monomial(q2 = 35.000000)}
					std::make_pair( 2226.4290588318418 , eps ),		// Im{P_GRV2026_monomial(q2 = 35.000000)}
					std::make_pair( 0.11 , eps ),		// Re{GRV2026_monomial_coeff_0}
					std::make_pair( -0.12 , eps ),		// Im{GRV2026_monomial_coeff_0}
					std::make_pair( 19.102252719737823 , eps ),		// Re{GRV2026_monomial_coeff_1}
					std::make_pair( 23.175962048833242 , eps ),		// Im{GRV2026_monomial_coeff_1}
					std::make_pair( -248.01713371245231 , eps ),		// Re{GRV2026_monomial_coeff_2}
					std::make_pair( -332.6185184216551 , eps ),		// Im{GRV2026_monomial_coeff_2}
					std::make_pair( 762.74820603631656 , eps ),		// Re{GRV2026_monomial_coeff_3}
					std::make_pair( 1251.8758951510599 , eps ),		// Im{GRV2026_monomial_coeff_3}
					std::make_pair( 819.61308726423977 , eps ),		// Re{GRV2026_monomial_coeff_4}
					std::make_pair( 798.42128372400714 , eps ),		// Im{GRV2026_monomial_coeff_4}
					std::make_pair( -4602.3254282863663 , eps ),		// Re{GRV2026_monomial_coeff_5}
					std::make_pair( -7826.1304798422616 , eps ),		// Im{GRV2026_monomial_coeff_5}
					std::make_pair( -0.00065870238051917698 , eps ),		// Re{H_+(q2 = -7.000000)}
					std::make_pair( 0.0013174047610383542 , eps ),		// Im{H_+(q2 = -7.000000)}
					std::make_pair( -0.00090315553726101931 , eps ),		// Re{H_+(q2 = -1.000000)}
					std::make_pair( 0.0010321777568697368 , eps ),		// Im{H_+(q2 = -1.000000)}
					std::make_pair( -0.12900744541590597 , eps ),		// Re{H_+(q2 = 5.000000)}
					std::make_pair( -0.19281787030543346 , eps ),		// Im{H_+(q2 = 5.000000)}
					std::make_pair( -1051.1156140136038 , eps ),		// Re{H_+(q2 = 16.000000)}
					std::make_pair( -3199.9373958569545 , eps ),		// Im{H_+(q2 = 16.000000)}
					std::make_pair( -27822.143563312977 , eps ),		// Re{H_+(q2 = 25.000000)}
					std::make_pair( 766.087099346017 , eps ),		// Im{H_+(q2 = 25.000000)}
					std::make_pair( 26060.480490101432 , eps ),		// Re{H_+(q2 = 35.000000)}
					std::make_pair( 7793.096993095055 , eps ),		// Im{H_+(q2 = 35.000000)}
					std::make_pair( 82429717.635367349 , eps ),		// GRV2026 weak bound
					std::make_pair( 86060987.502345353 , eps ),		// GRV2026 strong bound

                };
                TEST_CHECK_DIAGNOSTICS(diagnostics2, reference2);
            }
        }
} nonlocal_formfactor_grv2026_test;
#include <lib.hpp>


int main()
{
    function1D EvsP_B0, HvsP_B0;
    function2D E, P, M;
    std::cout << "Creating interpolators...";
    try {

        auto fileE = ".//eos//EvsH_B.dat";
        auto fileP = ".//eos//PvsH_B.dat";
        auto fileM = ".//eos//MvsH_B.dat";

        auto fileE0 = ".//eos//EvsP_B0.dat";
        auto fileH0 = ".//eos//HvsP_B0.dat";

        E = EoS::get_interpolator2D(fileE);
        P = EoS::get_interpolator2D(fileP);
        M = EoS::get_interpolator2D(fileM);

        EvsP_B0 = EoS::get_interpolator1D(fileE0);
        HvsP_B0 = EoS::get_interpolator1D(fileH0);

    }
    catch (const char * exp)
    {
        std::cout << "An error happened: " << exp << std::endl;
        return 0;
    }
    std::cout << "done." << std::endl;


    

    constexpr size_t N = 501;
    //constexpr double f0 = 4e-42;    
    std::ofstream file;
    file.open("sequence.dat");


    for ( double f0 = 0; f0 <= 2e-42; f0 += 0.5e-42 )
    {
        for (double Hc = 0.00001; Hc <= 0.024; Hc += 0.0005)
        {
            double Ec = E(Hc,0);
            double Pc = P(Hc,0);
            auto tov_guess = solve_tov(Pc,EvsP_B0,HvsP_B0, 100000*km, 0.01*km);
            auto solver =  PDEsolver(tov_guess,N,f0, E, P, M);

            try {
                solver.solve(Hc); 
                double cmf = solver.get_central_mf();
                double pmf = solver.get_polar_mf();
                double mass = solver.get_mass()/MS;
                double re = solver.get_equatorial_radius()/km;
                double rp = solver.get_polar_radius()/km;
                double pratio = solver.get_pressure_ratio();
                double Q = solver.get_quadrupolar_mass_moment()/MS/km/km;
                double me = solver.get_magnetic_energy()/MS;
                double Ec = E(Hc,std::abs(cmf));

                file << Ec << "\t" <<  Hc  <<  "\t"  << mass << "\t" << re << "\t" << rp << "\t" << Q << "\t" << me << "\t" << pratio <<
                    "\t" << cmf <<  "\t" << pmf  <<  "\t"  << f0 <<std::endl;
            }
            catch (...)
            {
                //std::cout << exp << std::endl;
                file.close();
            }
        }
    }
    
    
}
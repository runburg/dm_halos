////////////////////////////////////////////////////////////
// MADHAT (Model-Agnostic Dark Halo Analysis Tool )
//
// Inputs: takes in dwarf set, J factors, DM mass, and Inegrated Energy Spectrum
// Outputs: Bound on number of photons from dark matter for a confidence level, Phipp bound along, and DM cross section
// Author: Stephen Hill <hillsj@hawaii.edu>
// Date Last Modified: 07-22-2019
// Version Number: 0.4.0
////////////////////////////////////////////////////////////
#include <iostream>

#include <fstream>

#include <sstream>

#include <time.h>

#include <string>

#include <iomanip>

#include <boost/multiprecision/cpp_dec_float.hpp>

#include <boost/algorithm/string.hpp>

#include <boost/lexical_cast.hpp>

#include <boost/format.hpp>

using boost::multiprecision::cpp_dec_float;
typedef boost::multiprecision::number < cpp_dec_float < 50 >> madfloat;

using namespace std;
using namespace boost;
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
/////////////////////////////////////////////////////////////////////////////////////////////////
//All flags are false by defult to run the defult setup described in the README
bool fileout = true; //This flag changes the output for the first two output options to a file in output from outputing to the terminal
bool printall = true; //This flag changes the output so it displays beta as Nbound is incremented from outputing just the exact bound.
int gal = 47; //The number of dwarf galaxies in the data files
/////////////////////////////////////////////////////////////////////////////////////////////////
//dist() Finds the Probability of TOT DM
float dist(int Ndm, float Nbound) {
    float P2; //The float that will end out being the probability
    P2 = exp(Ndm * log(Nbound) - Nbound - lgamma(Ndm + 1));
    return P2; //Return the value
}
/////////////////////////////////////////////////////////////////////////////////////////////////
//READING DATA FROM DAT FILES
// Read Set
float process(float Beta, std::ifstream & INPUT, float mass, float energy, int Juse = 1, int Jerror = 1, int intype = 1) {
    float crosscons = 8 * M_PI * mass * mass / energy; //Creates a constant that is used to calculate the cross section
    int dwarf_count = 0;
    int dcol = 0;
    string line;
    string item;
    int header = 0;
    string skip("#");
    INPUT.seekg(0, ios::beg);
    while (getline(INPUT, line)) {
        if (contains(line, skip)) {
            header++;
        } else {
            break;
        }
    }
    INPUT.seekg(0, ios::beg);
    for (int i = 0; i < header + 1; ++i) {
        getline(INPUT, line);
    }
    while (getline(INPUT, line)) {
        ++dwarf_count;
        if (dwarf_count == 1) {
            stringstream ss(line);
            while (ss >> item)
                dcol++;
        }
    }
    INPUT.clear();
    INPUT.seekg(0, ios::beg);
    for (int i = 0; i < header + 1; ++i) {
        getline(INPUT, line);
    }
    float dwarf_data[dwarf_count][4];

    float check;
    for (int i = 0; i != dwarf_count; ++i) {
        dwarf_data[i][0] = 0;
        dwarf_data[i][1] = 0;
        dwarf_data[i][2] = 0;
        dwarf_data[i][3] = 0;
        getline(INPUT, line);
        istringstream read(line);
        while (read >> check) {
            read.seekg(0, ios::beg);
            read >> dwarf_data[i][0];
            read >> dwarf_data[i][1];
            read >> dwarf_data[i][2];
            read >> dwarf_data[i][3];
            break;
        }
    }
    INPUT.close();
    // Got Set_0
    // Get Obseravation data
    ifstream NOBS;
    NOBS.open("Data/NOBS.dat");
    if (!NOBS) {
        throw std::system_error(errno, std::system_category(), "failed to backend dwarf data");
    }
    header = 0;
    NOBS.seekg(0, ios::beg);
    while (getline(NOBS, line)) {
        if (contains(line, skip)) {
            header++;
        } else {
            break;
        }
    }
    NOBS.seekg(0, ios::beg);
    for (int i = 0; i < header; ++i) {
        getline(NOBS, line);
    }

    float obs_data[gal][3];
    for (int i = 0; i != gal; ++i) {
        for (int j = 0; j != 3; ++j) {
            NOBS >> obs_data[i][j];
        }
    }
    NOBS.close();
    // Got Observation Data
    // Get pmf
    ifstream PMF;
    PMF.open("Data/pmf.dat");
    if (!PMF) {
        throw std::system_error(errno, std::system_category(), "failed to backend fermi data");
    }
    header = 0;
    PMF.seekg(0, ios::beg);
    while (getline(PMF, line)) {
        if (contains(line, skip)) {
            header++;
        } else {
            break;
        }
    }
    PMF.seekg(0, ios::beg);
    for (int i = 0; i < header; ++i) {
        getline(PMF, line);
    }
    vector < vector < float >> pmf(20000, vector < float > (gal + 1));
    for (auto & i: pmf)
        std::fill(i.begin(), i.end(), 0);
    for (int i = 0; i != 861; ++i) {
        for (int j = 0; j != gal + 1; ++j) {
            PMF >> pmf[i][j];
        }
    }
    PMF.close();

    // Got pmf
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //Classifiying Dwarf Information
    int dwarf_list[dwarf_count]; //A list of dwarfs used by number
    for (int i = 0; i < dwarf_count; i++) {
        if (floorf(dwarf_data[i][0]) == dwarf_data[i][0] && dwarf_data[i][0] != 0 && dwarf_data[i][0] <= gal) //Checking if all dwarfs listed are grabable 
        {
            dwarf_list[i] = dwarf_data[i][0]; //importing the dwarf list for lookup
        } else {
            cout << "Not all dwarfs listed in set were findable please check integers against the list in the README" << endl;
            return 0;
        }
    }
    int Nobs = 0; //Number of observered over the set of dwarves
    for (int i = 0; i < dwarf_count; i++) {
        Nobs = Nobs + obs_data[dwarf_list[i] - 1][1]; //Summing over all observered to get the number
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //Calculating P1
    int Nbgd = 0;
    float P1[Nobs + 1];
    memset(P1, 0, sizeof(P1));
    if (dwarf_count == 1) { //For Number of dwarves=1 simply gathering them from the file
        for (int i = 0; i < (Nobs + 1); i++) {
            P1[i] = pmf[i][dwarf_list[0]];
        }
    }
    if (dwarf_count != 1) { //Incrementing up Ntot Bgd while summing over partitions
        float tempp = 0;
        float X[Nobs + 1];
        float I[Nobs + 1];

        for (int i = 0; i < Nobs + 1; i++) {
            I[i] = pmf[i][dwarf_list[0]];
        }
        for (int i = 0; i < Nobs + 1; i++) {
            X[i] = 0;
        }
        for (int k = 1; k < (dwarf_count - 1); k++) {
            for (int j = 0; j < Nobs + 2; j++) {
                for (int i = 0; i < j + 1; i++) {
                    tempp = pmf[i][dwarf_list[k]] * I[j - i];
                    X[j] = X[j] + tempp;
                }
            }
            for (int j = 0; j < Nobs + 1; j++) {
                I[j] = X[j];
            }
            for (int j = 0; j < Nobs + 1; j++) {
                X[j] = 0;
            }
        }
        for (int n = 0; n < (Nobs + 1); n++) {
            for (int k = 0; k < n + 1; k++) {
                X[k] = pmf[n - k][dwarf_list[dwarf_count - 1]] * I[k];
            }
            for (int k = 0; k < n + 1; k++) {
                P1[n] = P1[n] + X[k];
            }
            Nbgd++;
        }
    }
    Nbgd = 0;
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //Prepping J factors for calculation
    madfloat JAT = 0; //Variable for J factor times Aeff*Tobs
    madfloat JP = 0; //Variable for J factor times Aeff*Tobs + error
    madfloat JM = 0; //Variable for J factor times Aeff*Tobs -error
    madfloat SJAT = 0; //sum of J factors times Aeff*Tobs
    madfloat SJATP = 0; //+error sum of J factors times Aeff*Tobs
    madfloat SJATM = 0; //-error sum of J factors times Aeff*Tobs
    madfloat J = 0; //Jfactor
    if (Juse == 1) {
        for (int i = 0; i < dwarf_count; ++i) {
            J = pow(10, dwarf_data[i][1]); //Undoing the log base 10
            JAT = J * obs_data[dwarf_list[i] - 1][2]; //Multiplying by Aeff*Tobs
            SJAT = SJAT + JAT; //summing them up in the loop
            if (Jerror == 1) {
                JP = pow(10, dwarf_data[i][1]+dwarf_data[i][2]); //Taking the +error in log10(J) and undoing the log base 10
                JP = JP * obs_data[dwarf_list[i] - 1][2]; //Then muliplying by Aeff*Tobs
                SJATP = SJATP + JP; //summing those errors up
                // Then repeating the same process with M for minus
                JM = pow(10, dwarf_data[i][1]-dwarf_data[i][3]); //Taking the +error in log10(J) and undoing the log base 10
                JM = JM * obs_data[dwarf_list[i] - 1][2]; //Then muliplying by Aeff*Tobs
                SJATM = SJATM + JM; //summing those errors up
            }
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //Setting the headers for the output files
    if (Juse == 1 && Jerror == 1 && crosscons != 0 && mass == 0) //Has J factors Has J factor errors Has DM mass
    {
        printf("Nbound      Beta       Phi(cm^3 s^-1 GeV^-2)   +dPhi             -dPhi            cs(cm^3 s^-1)        +dcs             -dcs\n");
    }
    if (Juse == 1 && Jerror == 0 && crosscons != 0 && mass == 0) //Has J factors Has J factor errors Has DM mass
    {
        printf("Nbound      Beta       Phi(cm^3 s^-1 GeV^-2)   cs(cm^3 s^-1)\n");
    }
    if (Juse == 1 && Jerror == 1 && crosscons == 0 && mass == 0) //Has J factors Has J factor errors No DM mass
    {
        printf("Nbound      Beta       Phi(cm^3 s^-1 GeV^-2)   +dPhi             -dPhi\n");
    }
    if (Juse == 1 && Jerror == 0 && crosscons == 0 && mass == 0) //Has J factors No J factor errors No DM mass
    {
        printf("Nbound      Beta       Phi(cm^3 s^-1 GeV^-2)\n");
    }
    if (Juse == 0 && mass == 0) //No J factors No J factor errors No DM mass
    {
        printf("Nbound      Beta\n");
    }
    if (Juse == 1 && Jerror == 1 && crosscons != 0 && mass != 0) //Has J factors Has J factor errors Has DM mass
    {
        printf("Mass(Gev)   Energy       Nbound      Beta       Phi(cm^3 s^-1 GeV^-2)   +dPhi             -dPhi            cs(cm^3 s^-1)        +dcs             -dcs\n");
    }
    if (Juse == 1 && Jerror == 0 && crosscons != 0 && mass != 0) //Has J factors Has J factor errors Has DM mass
    {
        printf("Mass(Gev)   Energy       Nbound      Beta       Phi(cm^3 s^-1 GeV^-2)   cs(cm^3 s^-1)\n");
    }
    if (Juse == 1 && Jerror == 1 && crosscons == 0 && mass != 0) //Has J factors Has J factor errors No DM mass
    {
        printf("Mass(Gev)   Energy       Nbound      Beta       Phi(cm^3 s^-1 GeV^-2)   +dPhi             -dPhi\n");
    }
    if (Juse == 1 && Jerror == 0 && crosscons == 0 && mass != 0) //Has J factors No J factor errors No DM mass
    {
        printf("Mass(Gev)   Energy       Nbound      Beta       Phi(cm^3 s^-1 GeV^-2)\n");
    }
    if (Juse == 0 && mass != 0) //No J factors No J factor errors No DM mass
    {
        printf("Mass(Gev)      Energy       Nbound      Beta\n");
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //final summation loop
    float Pf = 0; //The complement of Beta
    float B = 0; //Beta  
    float P3;
    float PP;
    float Nbound = 1;
    for (Nbound = 1; Nbound < (Nobs + 1) && B < Beta; Nbound++) {
        int Ndm = 0;
        int Nbgd = Nobs; //Setting Ndm=Nobs
        Pf = 0; //Zeriong the placeholder for the sum
        float exp; //A variable to store the exponential
        for (int i = 0; i < (Nobs + 1); i++) //Summing over Nbgd+Ndm<NObs
        {
            P3 = 1;
            PP = 0;
            if (P1[Nbgd] > 0) {
                for (Ndm = 0; Ndm < (i + 1) && P3 > (0) || PP == 0; Ndm++) {
                    P3 = dist(Ndm, Nbound);
                    Pf = Pf + (P1[Nbgd] * P3); //Summing over the probability
                    PP = PP + P3;
                }
            }
            Nbgd--;
        }
        B = 1 - Pf; //Taking the complment of the sum
        if (printall == true && intype != 3) {
            madfloat PHI = Nbound; //A variable for PHI
            madfloat PHIP = Nbound; //A variable for PHI + error 
            madfloat PHIM = Nbound; //A variable for PHI - error
            madfloat CS; //A variable for CS
            madfloat CSP; //A variable for CS + error 
            madfloat CSM; //A variable for CS - error
            if (Juse == 1) {
                PHI = PHI / SJAT; //Phi=Nbound/JAT calculated earlier
            }
            if (Jerror == 1) {
                PHIP = PHIP / SJATP; //Phi=Nbound/JAT calculated earlier
                PHIM = PHIM / SJATM; //Phi=Nbound/JAT calculated earlier
                CS = PHI * crosscons;
                CSP = PHIP * crosscons;
                CSM = PHIM * crosscons;
            }
            if (mass != 0) {
                printf("%.4f     %.4f      ", mass, energy); //Print Nbound and Beta
            }
            printf("%.4f     %.4f      ", Nbound, B); //Print Nbound and Beta
            if (Juse == 1) {
                cout << PHI; //Print Phi if J factors
            }
            if (Jerror == 1) {
                printf("            ");
                cout << PHIM-PHI; //Print Phi error if J factors have error
                printf("        ");
                cout << -(PHIP-PHI);
            }
            if (crosscons != 0) {
                if (Juse == 1) {
                    printf("            ");
                    cout << CS; //Print Phi if J factors
                }
                if (Jerror == 1) {
                    printf("        ");
                    cout << CSM-CS; //Print Phi error if J factors have error
                    printf("        ");
                    cout << -(CSP-CS);
                }
            }
            printf("\n");
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //Final Calculation
    Nbound = Nbound - 2; //Incrementing back one to find the bound
    B = 0;
    //Loop counter that does the same as final summation but goes down decimals in percision
    for (int cut = 1; cut < 10; cut++) //cut is the decimal point
    {
        if (cut > 1) {
            Nbound = Nbound + pow(10, -cut - 1);
        }
        for (int acc = 0; acc < 9 && B < Beta; acc++) { //loops over .1,.2,.3...to .9
            int Ndm = 0;
            int Nbgd = Nobs; //Setting Ndm=Nobs
            Pf = 0; //Zeriong the placeholder for the sum
            float exp; //A variable to store the exponential
            for (int i = 0; i < (Nobs + 1); i++) //Summing over Nbgd+Ndm<NObs
            {
                P3 = 1;
                PP = 0;
                if (P1[Nbgd] > 0) {
                    for (Ndm = 0; Ndm < (i + 1) && P3 > (0) || PP == 0; Ndm++) {
                        P3 = dist(Ndm, Nbound);
                        Pf = Pf + (P1[Nbgd] * P3); //Summing over the probability
                        PP = PP + P3;
                    }
                }
                Nbgd--;
            }
            B = 1 - Pf; //Taking the complment of the sum
            Nbound = Nbound + pow(10, -cut);
        }
    }
    int Ndm = 0;
    Nbgd = Nobs; //Setting Ndm=Nobs
    Pf = 0; //Zeriong the placeholder for the sum
    float exp; //A variable to store the exponential
    for (int i = 0; i < (Nobs + 1); i++) //Summing over Nbgd+Ndm<NObs
    {
        P3 = 1;
        PP = 0;
        if (P1[Nbgd] > 0) {
            for (Ndm = 0; Ndm < (i + 1) && P3 > (0) || PP == 0; Ndm++) {
                P3 = dist(Ndm, Nbound);
                Pf = Pf + (P1[Nbgd] * P3); //Summing over the probability
                PP = PP + P3;
            }
        }
        Nbgd--;
    }
    B = 1 - Pf; //Taking the complment of the sum
    if (printall == false || intype == 3) {
        madfloat PHI = Nbound; //A variable for PHI
        madfloat PHIP = Nbound; //A variable for PHI + error 
        madfloat PHIM = Nbound; //A variable for PHI - error
        madfloat CS; //A variable for CS
        madfloat CSP; //A variable for CS + error 
        madfloat CSM; //A variable for CS - error
        if (Juse == 1) {
            PHI = PHI / SJAT; //Phi=Nbound/JAT calculated earlier
        }
        if (Jerror == 1) {
            PHIP = PHIP / SJATP; //Phi=Nbound/JAT calculated earlier
            PHIM = PHIM / SJATM; //Phi=Nbound/JAT calculated earlier
            CS = PHI * crosscons;
            CSP = PHIP * crosscons;
            CSM = PHIM * crosscons;
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////
        //Final Output
        if (mass != 0) {
            printf("%.4f     %.4f      ", mass, energy); //Print Nbound and Beta
        }
        printf("%.4f     %.4f      ", Nbound, B); //Print Nbound and Beta
        if (Juse == 1) {
            cout << PHI; //Print Phi if J factors
        }
        if (Jerror == 1) {
            printf("            ");
            cout << PHIM-PHI; //Print Phi error if J factors have error
            printf("        ");
            cout << -(PHIP-PHI);
        }
        if (crosscons != 0) {
            if (Juse == 1) {
                printf("        ");
                cout << CS; //Print Phi if J factors
            }
            if (Jerror == 1) {
                printf("        ");
                cout << CSM-CS; //Print Phi error if J factors have error
                printf("        ");
                cout << -(CSP-CS);
            }
        }
    }
    printf("\n");
    Beta = B;
    return Nbound;
}
/////////////////////////////////////////////////////////////////////////////////////////////////
// This is a duplicate of the main process if Nbound is already found which is the bulk of the calculation it exist for looping outputs
int output(float Nbound, float B, std::ifstream & INPUT, float mass, float energy, int Juse = 1, int Jerror = 1, int intype = 1) {
    float crosscons = 8 * M_PI * mass * mass / energy; //Creates a constant that is used to calculate the cross section   
    int dwarf_count = 0;
    int dcol = 0;
    string line;
    string item;
    int header = 0;
    string skip("#");
    INPUT.seekg(0, ios::beg);
    while (getline(INPUT, line)) {
        if (contains(line, skip)) {
            header++;
        } else {
            break;
        }
    }
    INPUT.seekg(0, ios::beg);
    for (int i = 0; i < header + 1; ++i) {
        getline(INPUT, line);
    }
    while (getline(INPUT, line)) {
        ++dwarf_count;
        if (dwarf_count == 1) {
            stringstream ss(line);
            while (ss >> item)
                dcol++;
        }
    }
    INPUT.clear();
    INPUT.seekg(0, ios::beg);
    for (int i = 0; i < header + 1; ++i) {
        getline(INPUT, line);
    }
    float dwarf_data[dwarf_count][4];

    float check;
    for (int i = 0; i != dwarf_count; ++i) {
        dwarf_data[i][0] = 0;
        dwarf_data[i][1] = 0;
        dwarf_data[i][2] = 0;
        dwarf_data[i][3] = 0;
        getline(INPUT, line);
        istringstream read(line);
        while (read >> check) {
            read.seekg(0, ios::beg);
            read >> dwarf_data[i][0];
            read >> dwarf_data[i][1];
            read >> dwarf_data[i][2];
            read >> dwarf_data[i][3];
            break;
        }
    }
    INPUT.close();
    // Got Set_0
    // Get Obseravation data
    ifstream NOBS;
    NOBS.open("Data/NOBS.dat");
    if (!NOBS) {
        throw std::system_error(errno, std::system_category(), "failed to backend dwarf data");
    }
    header = 0;
    NOBS.seekg(0, ios::beg);
    while (getline(NOBS, line)) {
        if (contains(line, skip)) {
            header++;
        } else {
            break;
        }
    }
    NOBS.seekg(0, ios::beg);
    for (int i = 0; i < header; ++i) {
        getline(NOBS, line);
    }
    float obs_data[gal][3];
    memset(obs_data, 0, sizeof(obs_data));
    for (int i = 0; i != gal; ++i) {
        for (int j = 0; j != 3; ++j) {
            NOBS >> obs_data[i][j];
        }
    }
    NOBS.close();
    // Got Observation Data
    // Get pmf
    ifstream PMF;
    PMF.open("Data/pmf.dat");
    if (!PMF) {
        throw std::system_error(errno, std::system_category(), "failed to backend fermi data");
    }
    header = 0;
    PMF.seekg(0, ios::beg);
    while (getline(PMF, line)) {
        if (contains(line, skip)) {
            header++;
        } else {
            break;
        }
    }
    PMF.seekg(0, ios::beg);
    for (int i = 0; i < header; ++i) {
        getline(PMF, line);
    }
    vector < vector < float >> pmf(20000, vector < float > (gal + 1));
    for (auto & i: pmf)
        std::fill(i.begin(), i.end(), 0);
    for (int i = 0; i != 861; ++i) {
        for (int j = 0; j != gal + 1; ++j) {
            PMF >> pmf[i][j];
        }
    }
    PMF.close();
    // Got pmf
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //Classifiying Dwarf Information
    int dwarf_list[dwarf_count]; //A list of dwarfs used by number
    for (int i = 0; i < dwarf_count; i++) {
        dwarf_list[i] = dwarf_data[i][0];
    }
    int Nobs = 0; //Number of observered over the set of dwarves
    for (int i = 0; i < dwarf_count; i++) {
        Nobs = Nobs + obs_data[dwarf_list[i] - 1][1]; //Summing over all observered to get the number
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //Calculating P1
    int Nbgd = 0;
    float P1[Nobs + 1];
    memset(P1, 0, sizeof(P1));
    if (dwarf_count == 1) { //For Number of dwarves=1 simply gathering them from the file
        for (int i = 0; i < (Nobs + 1); i++) {
            P1[i] = pmf[i][dwarf_list[0]];
        }
    }

    if (dwarf_count != 1) { //Incrementing up Ntot Bgd while summing over partitions
        float tempp = 0;
        float X[Nobs + 1];
        float I[Nobs + 1];

        for (int i = 0; i < Nobs + 1; i++) {
            I[i] = pmf[i][dwarf_list[0]];
        }
        for (int i = 0; i < Nobs + 1; i++) {
            X[i] = 0;
        }
        for (int k = 1; k < (dwarf_count - 1); k++) {
            for (int j = 0; j < Nobs + 2; j++) {
                for (int i = 0; i < j + 1; i++) {
                    tempp = pmf[i][dwarf_list[k]] * I[j - i];
                    X[j] = X[j] + tempp;
                }
            }
            for (int j = 0; j < Nobs + 1; j++) {
                I[j] = X[j];
            }
            for (int j = 0; j < Nobs + 1; j++) {
                X[j] = 0;
            }
        }
        for (int n = 0; n < (Nobs + 1); n++) {
            for (int k = 0; k < n + 1; k++) {
                X[k] = pmf[n - k][dwarf_list[dwarf_count - 1]] * I[k];
            }
            for (int k = 0; k < n + 1; k++) {
                P1[n] = P1[n] + X[k];
            }
            Nbgd++;
        }
    }
    Nbgd = 0;
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //Prepping J factors for calculation
    madfloat JAT = 0; //Variable for J factor times Aeff*Tobs
    madfloat JP = 0; //Variable for J factor times Aeff*Tobs + error
    madfloat JM = 0; //Variable for J factor times Aeff*Tobs -error
    madfloat SJAT = 0; //sum of J factors times Aeff*Tobs
    madfloat SJATP = 0; //+error sum of J factors times Aeff*Tobs
    madfloat SJATM = 0; //-error sum of J factors times Aeff*Tobs
    madfloat J = 0; //Jfactor
    if (Juse == 1) {
        for (int i = 0; i < dwarf_count; ++i) {
            J = pow(10, dwarf_data[i][1]); //Undoing the log base 10
            JAT = J * obs_data[dwarf_list[i] - 1][2]; //Multiplying by Aeff*Tobs
            SJAT = SJAT + JAT; //summing them up in the loop
            if (Jerror == 1) {
                JP = pow(10, dwarf_data[i][1]+dwarf_data[i][2]); //Undoing the log base 10
                JP = JP * obs_data[dwarf_list[i] - 1][2]; //Then muliplying by Aeff*Tobs
                SJATP = SJATP + JP; //summing those errors up
                // Then repeating the same process with M for minus
                JM = pow(10, dwarf_data[i][1]-dwarf_data[i][3]); //Taking the -error in log10(J) and undoing the log base 10
                JM = JM * obs_data[dwarf_list[i] - 1][2]; //Then muliplying by Aeff*Tobs
                SJATM = SJATM + JM; //summing those errors up
            }
        }
    }
    if (printall == false || intype == 3) {
        madfloat PHI = Nbound; //A variable for PHI
        madfloat PHIP= Nbound; //A variable for PHI + error 
        madfloat PHIM= Nbound; //A variable for PHI - error
        madfloat CS; //A variable for CS
        madfloat CSP; //A variable for CS + error 
        madfloat CSM; //A variable for CS - error
        if (Juse == 1) {
            PHI = PHI / SJAT; //Phi=Nbound/JAT calculated earlier
        }
        if (Jerror == 1) {
            PHIP = PHIP / SJATP; //Phi=Nbound/JAT calculated earlier
            PHIM = PHIM / SJATM; //Phi=Nbound/JAT calculated earlier
            CS = PHI * crosscons;
            CSP = PHIP * crosscons;
            CSM = PHIM * crosscons;
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////
        //Final Output
        if (mass != 0) {
            printf("%.4f     %.4f      ", mass, energy); //Print Nbound and Beta
        }
        printf("%.4f     %.4f      ", Nbound, B); //Print Nbound and Beta
        if (Juse == 1) {
            cout << PHI; //Print Phi if J factors
        }
        if (Jerror == 1) {
            printf("            ");
            cout << PHIM-PHI; //Print Phi error if J factors have error
            printf("        ");
            cout << -(PHIP-PHI);
        }
        if (crosscons != 0) {
            if (Juse == 1) {
                printf("        ");
                cout << CS; //Print Phi if J factors
            }
            if (Jerror == 1) {
                printf("        ");
                cout << CSM-CS; //Print Phi error if J factors have error
                printf("        ");
                cout << -(CSP-CS);
            }
        }
    }
    printf("\n");
    return 1;
}
/////////////////////////////////////////////////////////////////////////////////////////////////

int input() {
    //Telling it about your input file
    int Juse = 1; //Variable to know if there are J factors
    int Jerror = 1; //Variable to know if there are J factors have error
    char set[10]; //Input for which set to use
    char setJ[10]; //Input to know if there are J factors
    char setJE[10]; //Input to know if J factors have error
    ifstream INPUT;

    printf("Which set would you like to use?\n");
    fgets(set, 50, stdin);

    if (set[0] == '0') {
        INPUT.open("Set_0.dat");
        printf("Did you input J factors? Y/N\n");
        fgets(setJ, 50, stdin);

        if (setJ[0] == 'N' || setJ[0] == 'n' || setJ[1] == 'N' || setJ[1] == 'n') {
            Juse = 0;
            Jerror = 0;
        }
        if (Juse == 1) {
            printf("Did the J factors have error? Y/N\n");
            fgets(setJE, 50, stdin);
            if (setJE[0] == 'N' || setJE[0] == 'n' || setJE[1] == 'N' || setJE[1] == 'n') {
                Jerror = 0;
            }
        }
    }
    //Choosing from the defult types
    if (set[0] == '1') {
        printf("Which subset a,b,c,or Total\n");
        fgets(setJ, 50, stdin);
        if (setJ[0] == 'T' || setJ[0] == 't') {
            INPUT.open("DwarfSets/Set_1.dat");
        }
        if (setJ[0] == 'a' || setJ[0] == 'A') {
            INPUT.open("DwarfSets/Set_1a.dat");
        }
        if (setJ[0] == 'b' || setJ[0] == 'B') {
            INPUT.open("DwarfSets/Set_1b.dat");
        }
        if (setJ[0] == 'c' || setJ[0] == 'C') {
            INPUT.open("DwarfSets/Set_1c.dat");
        }
        if (setJ[0] != 'T' && setJ[0] != 't' && setJ[0] != 'A' && setJ[0] != 'a' && setJ[0] != 'B' && setJ[0] != 'b' && setJ[0] != 'C' && setJ[0] != 'c') {
            cout << "not a valid Dwarf Set" << endl;
            return 0;
        }
    }
    if (set[0] == '2') {
        INPUT.open("DwarfSets/Set_2.dat");
    }
    if (set[0] == '3') {
        INPUT.open("DwarfSets/Set_3.dat");
    }
    if (set[0] == '4') {
        INPUT.open("DwarfSets/Set_4.dat");
    }
    if (set[0] == '5') {
        INPUT.open("DwarfSets/Set_5.dat");
    }
    if (set[0] != '5' && set[0] != '4' && set[0] != '3' && set[0] != '2' && set[0] != '1' && set[0] != '0') {
        cout << "Not A valid Dwarf Set" << endl;
        return 0;
    }
    printf("Specify a Confidence Level for the results.  Eg. 0.68, 0.90, 0.95, etc.\n");
    float Beta;
    cin >> Beta;
    if (!cin || Beta >= 1 || Beta < .001) // or if(cin.fail())
    {
        cout << "Not A valid Beta Value" << endl;
        return 0;
    }
    printf("Specify the dark matter mass in GeV, or 0 (zero) if not applicable.\n");
    float mass;
    cin >> mass;
    if (!cin || mass < 0) // or if(cin.fail())
    {
        cout << "Not A valid Mass Value" << endl;
        return 0;
    }
    float energy = 1;
    if (mass != 0) {
        printf("Specify the integrated photon spectrum between 1 GeV and 100 GeV.\n");
        float energy;
        cin >> energy;
        if (!cin || energy <= 0) // or if(cin.fail())
        {
            cout << "Not A valid Energy Value" << endl;
            return 0;
        }
    }
    float crosscons = 8 * M_PI * mass * mass / energy; //Creates a constant that is used to calculate the cross section
    if (fileout == true) {
        string outst = "Output/DM_"; //create output string
        outst = outst + boost::lexical_cast < std::string > (boost::format("%.4f") % Beta); //Add confidence
        time_t timer;
        outst = outst + boost::lexical_cast < std::string > (boost::format("%.4f") % time( & timer)); //Add confidence
        outst += ".out"; //Add .out
        cout << "File ";
        cout << outst;
        cout << " created " << endl;
        const char * out = outst.c_str(); //setting the output string to the file name
        stdout = freopen(out, "w", stdout); // Setting the the print to the output file
    }
    process(Beta, INPUT, mass, energy, Juse, Jerror, 1);
}
/////////////////////////////////////////////////////////////////////////////////////////////////
int input(char * dat, float Beta) {
    //Telling it about your input file
    int Juse = 1; //Variable to know if there are J factors
    int Jerror = 1; //Variable to know if there are J factors have error
    char set[10]; //Input for which set to use
    char setJ[10]; //Input to know if there are J factors
    char setJE[10]; //Input to know if J factors have error
    ifstream INPUT;
    INPUT.open(dat);
    if (!INPUT)
        throw std::system_error(errno, std::system_category(), "failed to open dwarf set");
    printf("Specify the dark matter mass in GeV, or 0 (zero) if not applicable.\n");
    float mass;
    cin >> mass;
    if (!cin || mass < 0) // or if(cin.fail())
    {
        cout << "Not A valid Mass Value" << endl;
        return 0;
    }
    float energy = 1;
    if (mass != 0) {
        printf("Specify the integrated photon spectrum between 1 GeV and 100 GeV.\n");
        float energy;
        cin >> energy;
        if (!cin || energy <= 0) // or if(cin.fail())
        {
            cout << "Not A valid Energy Value" << endl;
            return 0;
        }
    }
    float crosscons = 8 * M_PI * mass * mass / energy; //Creates a constant that is used to calculate the cross section
    if (fileout == true) {
        string outst = "Output/DM_"; //create output string
        string name;
        name = dat;
        boost::erase_all(name, "/");
        outst += name; //Add set name
        outst = outst.substr(0, outst.size() - 4); //remove set file extension
        outst += "_"; //Add underscore
        outst = outst + boost::lexical_cast < std::string > (boost::format("%.4f") % Beta); //Add confidence
        time_t timer;
        outst = outst + boost::lexical_cast < std::string > (boost::format("%.4f") % time( & timer)); //Add confidence
        outst += ".out"; //Add .out
        cout << "File ";
        cout << outst;
        cout << " created " << endl;
        const char * out = outst.c_str(); //setting the output string to the file name
        stdout = freopen(out, "w", stdout); // Setting the the print to the output file
    }
    process(Beta, INPUT, mass, energy, Juse, Jerror, 2);
}
/////////////////////////////////////////////////////////////////////////////////////////////////
int input(char * dat, float Beta, char * in ) {
    //Telling it about your input file
    //Reading the in file
    int line_count = 0;
    int dcol = 0;
    string line;
    string item;
    ifstream INPUT;
    ifstream MASS;
    MASS.open( in );
    if (!MASS)
        throw std::system_error(errno, std::system_category(), "failed to open mass and spectrum data");
    int header = 0;
    string skip("#");
    MASS.seekg(0, ios::beg);
    while (getline(MASS, line)) {
        if (contains(line, skip)) {
            header++;
        } else {
            break;
        }
    }
    MASS.seekg(0, ios::beg);
    for (int i = 0; i < header + 1; ++i) {
        getline(MASS, line);
    }
    while (getline(MASS, line)) {
        ++line_count;
        if (line_count == 1) {
            stringstream ss(line);
        }
    }
    MASS.clear();
    MASS.seekg(0, ios::beg);
    for (int i = 0; i < header + 1; ++i) {
        getline(MASS, line);
    }
    float mass_data[line_count][2];
    memset(mass_data, 0, sizeof(mass_data));
    float check;
    for (int i = 0; i != line_count; ++i) {
        mass_data[i][0] = 0;
        mass_data[i][1] = 0;
        getline(MASS, line);
        istringstream read(line);
        while (read >> check) {
            read.seekg(0, ios::beg);
            read >> mass_data[i][0];
            read >> mass_data[i][1];
            break;
        }
    }
    MASS.close();
    string name;
    int Juse = 1; //Variable to know if there are J factors
    int Jerror = 1; //Variable to know if there are J factors have error
    char set[10]; //Input for which set to use
    char setJ[10]; //Input to know if there are J factors
    char setJE[10]; //Input to know if J factors have error
    string outst = "Output/DM_"; //create output string
    name = dat;
    boost::erase_all(name, "/");
    outst += name; //Add set name
    outst = outst.substr(0, outst.size() - 4); //remove set file extension
    outst += "_"; //Add underscore
    outst = outst + boost::lexical_cast < std::string > (boost::format("%.4f") % Beta); //Add confidence
    outst += "_"; //Add underscore
    time_t timer;
    outst = outst + boost::lexical_cast < std::string > (boost::format("%.4f") % time( & timer)); //Add confidence
    outst += ".out"; //Add .out
    cout << "File ";
    cout << outst;
    cout << " created " << endl;
    const char * out = outst.c_str(); //setting the output string to the file name
    stdout = freopen(out, "w", stdout); // Setting the the print to the output file
    INPUT.open(dat);
    if (!INPUT)
        throw std::system_error(errno, std::system_category(), "failed to open dwarf set");
    float mass = mass_data[0][0];
    float energy = mass_data[0][1];
    float Nbound = process(Beta, INPUT, mass, energy, Juse, Jerror, 3); //Getting Nbound from the calculation
    for (int loop = 1; loop < line_count; ++loop) {
        INPUT.open(dat);
        float mass = mass_data[loop][0]; //Changing the mass
        float energy = mass_data[loop][1]; //Changing the energy
        output(Nbound, Beta, INPUT, mass, energy, Juse, Jerror, 3); //Running just the output calculations
    }
}
/////////////////////////////////////////////////////////////
int main(int argc, char * argv[]) {
    cout << "\n Model-Agnostic Dark Halo Analysis Tool  \n";
    if (argc == 1 || argc == 2) {
        input(); //Run The Main Function
    }
    if (argc == 3) {
        float Beta = atof(argv[2]);
        if (Beta > 0 && Beta < 1) {
            input(argv[1], Beta);
        } else {
            cout << "Invalid Beta Value" << endl;
        }
    }
    if (argc == 4) {
        float Beta = atof(argv[2]);
        if (Beta > 0 && Beta < 1) {
            input(argv[1], Beta, argv[3]);
        } else {
            cout << "Invalid Beta Value" << endl;
        }
    }
    return 0; //End The Program
}

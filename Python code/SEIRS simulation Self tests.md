```python
# This code uses the SEIRS simulation function to investigate the transmission of COVID-19 in the workplace if tests had
# the sensitivity as stated in Stohr et al. (2021).
```


```python
# Import all the necessary packages

from seirsplus.models import *
from seirsplus.networks import *
#from seirsplus.sim_loops import *
from seirsplus.utilities import *
import networkx
import matplotlib.pyplot as pyplot
import seaborn
import pandas as pd
import numpy as np
import os
import math
import sys
```


```python
# Set working directory to location where files are located
os.chdir("C:\\Users\\Quint\\Documents\\Documents\\Master Behavioral Data Science\\MasterThese\\data\\seirsplus_package_thesis\\build\\lib\\seirsplus")

# Import self adapted code of the seirsplus package
# Source code altered to implement 'none' PCR Frequency and with changed sensitivity of PCR tests according to Stohr et al. (2021).
from sim_loops_altered_self_test import *
```


```python
# Delete all dataframes
del values_df
del combined_df
```


```python
# Create empty dataframe to store values per condition
values_df = pd.DataFrame(columns=['Network_size','Bubble_size','PCR_frequency','R','Percentage_infected', 
                                  'Peak_Percentage_Hospitalized', "Percentage_Fatality"])

```


```python
# Create dataframe to store all values of every condition in
combined_df = pd.DataFrame()
```


```python
# Function to introduce 1 infection per 2 months for net size = 10
# Use same proportion for bigger network sizes
def average_introduction(Network_size):
    if Network_size == 20:
        global introduction
        introduction = 1/60
    elif Network_size == 50:
        introduction = 2.5/60
    elif Network_size == 100:
        introduction = 5/60
    elif Network_size == 250:
        introduction = 12.5/60
    elif Network_size == 500:
        introduction = 25/60
    elif Network_size == 1000:
        introduction = 50/60
```


```python
# Constructing simulation function
def simulation_covid(Network_size, Bubble_size, PCR_Frequency, R, Itterations):

    # Calculate introduction for every itteration
    average_introduction(Network_size)
        
    # Specify loop for iterations per condition
    for i in range(Itterations):
            
            # Print parameters
            print('Network size =', Network_size)
            print('Bubble size =', Bubble_size)
            print('PCR_frequency =', PCR_Frequency)
            print('R =', R)
        
        
            # Set basic parameters
            NUM_COHORTS              = Network_size//Bubble_size
            NUM_NODES_PER_COHORT     = Bubble_size
            NUM_TEAMS_PER_COHORT     = 1

            MEAN_INTRACOHORT_DEGREE  = round(0.10*Bubble_size)  
                                                                    
            PCT_CONTACTS_INTERCOHORT = 0.0

            N = NUM_NODES_PER_COHORT*NUM_COHORTS
            # Initial prevalence of the disease in the network
            INIT_EXPOSED = (.10*Network_size)            # Not realistic, but all networks have same % infected


            # Farz network
            # alpha = strength of common neighbor's effect on edge formation (tunes transitivity, clustering)
            # gamma = strength of degree similarity effect on edge formation (tunes assortativity)
            # beta = probability of edges formation within communities, rather than between (strength of community structure)
            # r = maximum number of communities each node can belong to
            # q = probability of a node belonging to the multiple communities
            # phi = constant added to all community sizes, higher number makes the communities more balanced in size, 1 results in power law community size distribution
            # epsilon = probability of noisy/random edges
            G_baseline, cohorts, teams = generate_workplace_contact_network(
                                             num_cohorts=NUM_COHORTS, num_nodes_per_cohort=NUM_NODES_PER_COHORT,
                                             num_teams_per_cohort=NUM_TEAMS_PER_COHORT,
                                             mean_intracohort_degree=MEAN_INTRACOHORT_DEGREE,
                                             pct_contacts_intercohort=PCT_CONTACTS_INTERCOHORT,
                                             farz_params={'alpha':5.0, 'gamma':5.0, 'beta':0.5, 'r':1, 'q':0.0, 'phi':10,
                                                          'b':0, 'epsilon':1e-6, 'directed': False, 'weighted': False})

                        # Plot (if you want)
            #network_info(G_baseline, "Baseline", plot=True)

            # Here we define quarantine contact network to be an empty one (no connections)
            G_quarantine = networkx.classes.function.create_empty_copy(G_baseline)

            # Here we generate distributions of values for each parameter, thus specifying a realistic hetereogeneous population
            # For the latent period
            latentPeriod_mean, latentPeriod_coeffvar = 2.5, 0.6
            SIGMA   = 1 / gamma_dist(latentPeriod_mean, latentPeriod_coeffvar, N)
            # For the pre-symptomatic period
            presymptomaticPeriod_mean, presymptomaticPeriod_coeffvar = 3.5, 0.6
            LAMDA   = 1 / gamma_dist(presymptomaticPeriod_mean, presymptomaticPeriod_coeffvar, N)
            # Make the plot where latent period, pre-symptomatic period & total incubation period are shown (if you want)
            #dist_info([1/LAMDA, 1/SIGMA, 1/LAMDA+1/SIGMA], ["latent period", "pre-symptomatic period", "total incubation period"],
             #         plot=True, colors=['gold', 'darkorange', 'black'], reverse_plot=True)


            # Here we generate distributions expected (a)symptomatic periods (time in symptomatic or asymptomatic state).
            symptomaticPeriod_mean, symptomaticPeriod_coeffvar = 6.5, 0.4
            GAMMA   = 1 / gamma_dist(symptomaticPeriod_mean, symptomaticPeriod_coeffvar, N)
            # The expected total infectious period for each individual is the sum of their expected pre-symptomatic and 
            # (a)symptomatic periods.
            # Specify infectious period
            infectiousPeriod = 1/LAMDA + 1/GAMMA
            # Make the plot (if you want)
            #dist_info([1/LAMDA, 1/GAMMA, 1/LAMDA+1/GAMMA], ["pre-symptomatic period", "(a)symptomatic period", "total infectious period"],
                      #plot=True, colors=['darkorange', 'crimson', 'black'], reverse_plot=True)



            # Generate distribution of expected onset to hospitalization periods
            onsetToHospitalizationPeriod_mean, onsetToHospitalizationPeriod_coeffvar = 8.0, 0.45
            ETA     = 1 / gamma_dist(onsetToHospitalizationPeriod_mean, onsetToHospitalizationPeriod_coeffvar, N)
            # Hospitalization discharge
            hospitalizationToDischargePeriod_mean, hospitalizationToDischargePeriod_coeffvar = 14.0, 0.45
            GAMMA_H = 1 / gamma_dist(hospitalizationToDischargePeriod_mean, hospitalizationToDischargePeriod_coeffvar, N)
            # Plot (if you want)
            #dist_info([1/ETA, 1/GAMMA_H, 1/ETA+1/GAMMA_H], ["onset-to-hospitalization period", "hospitalization-to-discharge period", "onset-to-discharge period"],
                      #plot=True, colors=['crimson', 'violet', 'black'], reverse_plot=True)

            # Generate distribution of hospitalization to death
            hospitalizationToDeathPeriod_mean, hospitalizationToDeathPeriod_coeffvar = 7.0, 0.45
            MU_H    = 1 / gamma_dist(hospitalizationToDeathPeriod_mean, hospitalizationToDeathPeriod_coeffvar, N)
            # Plot (if you want)
            #dist_info([1/ETA, 1/MU_H, 1/ETA+1/MU_H], ["onset-to-hospitalization period", "hospitalization-to-death period", "onset-to-death period"],
                    #  plot=True, colors=['crimson', 'darkgray', 'black'], reverse_plot=True)

            # Severity parameters
            # This percentage of case will progress from the pre-symptomatic state to the asymptomatic state, 
            # rather than to the symptomatic state.
            PCT_ASYMPTOMATIC = 0.308     
            # Hospitalization rate
            PCT_HOSPITALIZED = 0.043
            # Fatality rate for hospitalizations
            PCT_FATALITY = 0.05

            # Transmission parameters
            # Mean expected number of new infections by a single infectious individual
            R0_mean     = R
            # How much the mean varies, low is not much variance
            R0_coeffvar = 0.2
            # Expected number of new infections by a single infectious individual
            R0 = gamma_dist(R0_mean, R0_coeffvar, N)
            # Plot (if you want)
            #dist_info(R0, "Individual R0", bin_size=0.1, plot=True, colors='crimson')

            # Individuals get an Individual Tranmitability Parameter stored in BETA
            BETA = 1/infectiousPeriod * R0

            # 40% of interactions being with incidental or casual contacts outside their set of close contacts(Co-workers)
            P_GLOBALINTXN = 0.4

            # Testing, Tracing & Isolation intervention protocal measures
            INTERVENTION_START_PCT_INFECTED = 0/100
            
            # Expected number of new exogenous exposures per day
            AVERAGE_INTRODUCTIONS_PER_DAY   = introduction  #(0.00026470588*Network_size)       

            # How often to do testing (other than self-reporting symptomatics who can get tested any day)
            TESTING_CADENCE                 = PCR_Frequency 
            
            # Max daily test allotment defined as a percent of population size
            PCT_TESTED_PER_DAY              = 1.0         
            
            # Test false negative rate, will use FN rate that varies with disease time
            TEST_FALSENEG_RATE              = 'temporal'   
            
             # Max percent of daily test allotment to use on self-reporting symptomatics
            MAX_PCT_TESTS_FOR_SYMPTOMATICS  = 1.0         
            
            # Max percent of daily test allotment to use on contact traces
            MAX_PCT_TESTS_FOR_TRACES        = 0.0           
            
            # Magnitude of degree bias in random selections for testing, none here
            RANDOM_TESTING_DEGREE_BIAS      = 0             

             # Percentage of primary cases' contacts that are traced
            PCT_CONTACTS_TO_TRACE           = 0.0           
            
            # Number of cadence testing days between primary tests and tracing tests
            TRACING_LAG                     = 2             

            
            # Number of days between onset of symptoms and self-isolation of symptomatics
            ISOLATION_LAG_SYMPTOMATIC       = 2             
            
            # Test turn-around time (TAT): number of days between administration of test and isolation of positive cases
            ISOLATION_LAG_POSITIVE          = 1             
            
            # Number of days between a contact being traced and that contact self-isolating
            ISOLATION_LAG_CONTACT           = 0             

            # Intervention compliance parameters
            TESTING_COMPLIANCE_RATE_SYMPTOMATIC                  = 0.75
            TESTING_COMPLIANCE_RATE_TRACED                       = 0.0
            
            # Assume employee testing is mandatory, so 100% compliance
            TESTING_COMPLIANCE_RATE_RANDOM                       = 1.0  
            TRACING_COMPLIANCE_RATE                              = 0.0
            ISOLATION_COMPLIANCE_RATE_SYMPTOMATIC_INDIVIDUAL     = 0.0
            ISOLATION_COMPLIANCE_RATE_SYMPTOMATIC_GROUPMATE      = 0.0
            ISOLATION_COMPLIANCE_RATE_POSITIVE_INDIVIDUAL        = 1.0
            ISOLATION_COMPLIANCE_RATE_POSITIVE_GROUPMATE         = 0.0
            
            # 100% of those who have been in contact with an infected person, will go in isolation
            ISOLATION_COMPLIANCE_RATE_POSITIVE_CONTACT           = 1.0 
            ISOLATION_COMPLIANCE_RATE_POSITIVE_CONTACTGROUPMATE  = 0.0

            # Now we randomly assign true/false compliance to each individual according to the rates above. 
            # True will participate, False will not
            TESTING_COMPLIANCE_RANDOM                        = (numpy.random.rand(N) < TESTING_COMPLIANCE_RATE_RANDOM)
            TESTING_COMPLIANCE_TRACED                        = (numpy.random.rand(N) < TESTING_COMPLIANCE_RATE_TRACED)
            TESTING_COMPLIANCE_SYMPTOMATIC                   = (numpy.random.rand(N) < TESTING_COMPLIANCE_RATE_SYMPTOMATIC)

            TRACING_COMPLIANCE                               = (numpy.random.rand(N) < TRACING_COMPLIANCE_RATE)

            ISOLATION_COMPLIANCE_SYMPTOMATIC_INDIVIDUAL      = (numpy.random.rand(N) < ISOLATION_COMPLIANCE_RATE_SYMPTOMATIC_INDIVIDUAL)
            ISOLATION_COMPLIANCE_SYMPTOMATIC_GROUPMATE       = (numpy.random.rand(N) < ISOLATION_COMPLIANCE_RATE_SYMPTOMATIC_GROUPMATE)
            ISOLATION_COMPLIANCE_POSITIVE_INDIVIDUAL         = (numpy.random.rand(N) < ISOLATION_COMPLIANCE_RATE_POSITIVE_INDIVIDUAL)
            ISOLATION_COMPLIANCE_POSITIVE_GROUPMATE          = (numpy.random.rand(N) < ISOLATION_COMPLIANCE_RATE_POSITIVE_GROUPMATE)
            ISOLATION_COMPLIANCE_POSITIVE_CONTACT            = (numpy.random.rand(N) < ISOLATION_COMPLIANCE_RATE_POSITIVE_CONTACT)
            ISOLATION_COMPLIANCE_POSITIVE_CONTACTGROUPMATE   = (numpy.random.rand(N) < ISOLATION_COMPLIANCE_RATE_POSITIVE_CONTACTGROUPMATE)

            # Making the model
            model = ExtSEIRSNetworkModel(G=G_baseline, p=P_GLOBALINTXN,
                                          beta=BETA, sigma=SIGMA, lamda=LAMDA, gamma=GAMMA,
                                          gamma_asym=GAMMA, eta=ETA, gamma_H=GAMMA_H, mu_H=MU_H,
                                          a=PCT_ASYMPTOMATIC, h=PCT_HOSPITALIZED, f=PCT_FATALITY,
                                          G_Q=G_quarantine, isolation_time=14,
                                          initE=INIT_EXPOSED)
            
            # Max simulation time to 100 days
            # If you want this to be higher, change testing cadence & cadence_cycle_length in source code: sim_loops_altered
            # Make it slightly higher than T
            T = 100

            # Execute the TTI simulatiom
            run_tti_sim(model, T,
                        intervention_start_pct_infected=INTERVENTION_START_PCT_INFECTED, average_introductions_per_day=AVERAGE_INTRODUCTIONS_PER_DAY,
                        testing_cadence=TESTING_CADENCE, pct_tested_per_day=PCT_TESTED_PER_DAY, test_falseneg_rate=TEST_FALSENEG_RATE,
                        testing_compliance_symptomatic=TESTING_COMPLIANCE_SYMPTOMATIC, max_pct_tests_for_symptomatics=MAX_PCT_TESTS_FOR_SYMPTOMATICS,
                        testing_compliance_traced=TESTING_COMPLIANCE_TRACED, max_pct_tests_for_traces=MAX_PCT_TESTS_FOR_TRACES,
                        testing_compliance_random=TESTING_COMPLIANCE_RANDOM, random_testing_degree_bias=RANDOM_TESTING_DEGREE_BIAS,
                        tracing_compliance=TRACING_COMPLIANCE, pct_contacts_to_trace=PCT_CONTACTS_TO_TRACE, tracing_lag=TRACING_LAG,
                        isolation_compliance_symptomatic_individual=ISOLATION_COMPLIANCE_SYMPTOMATIC_INDIVIDUAL, isolation_compliance_symptomatic_groupmate=ISOLATION_COMPLIANCE_SYMPTOMATIC_GROUPMATE,
                        isolation_compliance_positive_individual=ISOLATION_COMPLIANCE_POSITIVE_INDIVIDUAL, isolation_compliance_positive_groupmate=ISOLATION_COMPLIANCE_POSITIVE_GROUPMATE,
                        isolation_compliance_positive_contact=ISOLATION_COMPLIANCE_POSITIVE_CONTACT, isolation_compliance_positive_contactgroupmate=ISOLATION_COMPLIANCE_POSITIVE_CONTACTGROUPMATE,
                        isolation_lag_symptomatic=ISOLATION_LAG_SYMPTOMATIC, isolation_lag_positive=ISOLATION_LAG_POSITIVE,
                        isolation_groups=list(teams.values()))
            
            
            # Test for calculating percentage
            percentages = ((model.total_num_infected()[-1]+model.total_num_recovered()[-1])/model.numNodes * 100) 
            # Total percent fatality: 
            fatality = (model.numF[-1]/model.numNodes * 100)
            # Peak  pct hospitalized: 
            hospitalized = (numpy.max(model.numH)/model.numNodes * 100)

            # Store values of itteration in temporary dataframe
            values_df.loc[i, ['Network_size']] = NUM_COHORTS*NUM_NODES_PER_COHORT
            values_df.loc[i, ['Bubble_size']] = Bubble_size
            values_df.loc[i, ['PCR_frequency']] = PCR_Frequency
            values_df.loc[i, ['R']] = R
            values_df.loc[i, ['Percentage_infected']] = percentages
            values_df.loc[i, ['Peak_Percentage_Hospitalized']] = hospitalized
            values_df.loc[i, ["Percentage_Fatality"]] = fatality
                        
            #Show results
            results_summary(model)
            
            # Print to indicate next iteration is starting
            print('-----------------------------------------------------')

            # Store temporary values in combined dataframe
            global combined_df
            combined_df = combined_df.append(values_df, True)


        

        


```


```python
# Loop to add the no bubble conditions
Network_size = [20, 50, 100, 250]
PCR_Frequency = ['workday', 'weekly', 'semiweekly', 'monthly', 'none']
R = [1, 1.5, 2, 2.5]

for i in range(len(Network_size)):
    for k in range(len(PCR_Frequency)):
        for l in range(len(R)):
            simulation_covid(Network_size[i], Network_size [i], PCR_Frequency[k], R[l], 100)
```


```python
# Add none for bubble size if bubble size = network size
combined_df.loc[(combined_df.Bubble_size == combined_df.Network_size),'Bubble_size'] = 'none'
```


```python
# Check if loop succeeded
combined_df
```


```python
# Change rest of the colnames for an ordered dataframe
combined_df = combined_df.rename(columns = {"Network_size":"Network size", "Bubble_size":"Bubble size",
                                           "PCR_frequency":"PCR frequency", "Percentage_infected":"Percentage infected",
                                            "Peak_Percentage_Hospitalized":"Peak percentage hospitalized",
                                            "Percentage_Fatality": "Percentage fatality"})

```


```python
# Calculate the mean and SD per condition for percentage infected
combined_df['Percentage infected'] = combined_df['Percentage infected'].astype(np.int64)

# Create a new dataframe for the means and SD's
final_df = combined_df.groupby(['Network size', 'Bubble size', 'PCR frequency', 'R']).agg({'Percentage infected': 
                                                                                              ['mean', 'std']})
```


```python
# Make final data frame
final_df.columns = ['mean', 'sd']
final_df = final_df.reset_index()
final_df
```


```python
# Save final data frame as Final data as excel
final_df.to_csv(r'C:\Users\Quint\Documents\Documents\Master Behavioral Data Science\MasterThese\Data\Full_dataset_selftest.csv')

# Save as combined data as combined data as excel
combined_df.to_csv(r'C:\Users\Quint\Documents\Documents\Master Behavioral Data Science\MasterThese\Data\Combined_dataset_selftest.csv')
```


```python

```

# Tropical Cyclones (hk25-TropCyc)

Tropical Cyclones (TCs) are intense organised convective systems that are responsible for nearly half of the worldwide disaster-related costs (CRED & UNDRR, 2020). Historically, they have been difficult to represent in coarse-resolution climate models because of their small size and their sensitivity to convective parametrisations. It has been shown that increasing resolution to 25km allows models to represent the number and distribution of cyclones correctly (Roberts et al. 2020), but the intensity remains largely underestimated, and the structure is not well represented (Bourdin et al. 2024). Baker et al. (2024) showed that increasing resolution up to 5km improves the realism of intensity, intensification rate, lifecycle and structure of TCs in NextGEMS simulations.

In this team, we investigated how TCs are represented in the new sets of simulations in terms of statistics, structure, lifecycle and link with the environment.

**Coordination**: Stella Bourdin (stella.bourdin@physics.ox.ac.uk), Alex Baker (alexander.baker@reading.ac.uk), Arthur Avenas (arthur.avenas@esa.int), Xu Chen (chenx@g.ecc.u-tokyo.ac.jp)

In this folder we leave for legacy scripts to track TCs in the simulations, as well as some of the data (in particular track datasets).

## Useful resources

* [TempestExtremes](https://github.com/ClimateGlobalChange/tempestextremes)
* [HuracanPy](https://huracanpy.readthedocs.io/en/latest/)

## Connectivity Files for TempestExtremes

With great help from Bryce Harrop, we are now able to run TempestExtremes on HealPix, using dedicated ConnectivityFiles. 
Files can be generated using instructions in `notebooks/1_Get_Connectivity_Files/`. Corresponding files for zoom 0-2 are provided via GitHub in `ConnectivityFiles` while files for zoom 3-10 can be found on JASMIN in `/home/users/sbourdin/WCRP_Hackathon/hk25-teams/hk25-TropCyc/ConnectivityFiles/` or in [this Google Drive folder](https://drive.google.com/drive/folders/1fNDDQA_G-yy05SP8J8pV2EIFJCs1bOtb?usp=sharing).

## TC tracking



3. Get familiar with the tracking process: Check-out the pre-processing, TE_intro and TE_full notebooks to understand how it works.

4. Run the tracking with the full_tracking notebook (or turn it into a script if necessary).

5. Check the tracks you obtained (An image will be generated at the end of the full_tracking script to get a first glance.

6. Upload the tracks to hk25-TropCyc/TC_tracks/ with one csv file per simulation, named with the simulation code. Update the Monitoring TC tracks once this is done

*Technical details*:
TCs are to be tracked in all relevant simulations (i.e. simulations covering TC-prone areas) using TempestExtremes and the algorithm described in Ullrich et al. 2021.
Output should include, a minima, track_id, time, lon, lat, maximum wind speed in a 2°GCD radius and minimum SLP.

Note: If some simulations don't have the necessary data for the tracking, an alternative 2D algorithm can be applied - please reach out if that's the case.

![image](tracking_process.png)
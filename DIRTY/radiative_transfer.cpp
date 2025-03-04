// ======================================================================
//   Procedure to do Monte Carlo radiative transfer for photons through
// dust.  This procedure assumes the dust grid is already constructed.
//
// 2004 Dec/KDG - written
// 2005 May/KDG - added output stuff
// 22 Oct 2013/KDG - added option to deal with forced_first_scatter getting zero
// optical depth
// ======================================================================
#include "radiative_transfer.h"
// #define DEBUG_RT
// #define OUTNUM -1
// #define OUTNUM 453
#define DEBUG_OUTRANGE

void radiative_transfer(geometry_struct &geometry, runinfo_struct &runinfo, output_struct &output, photon_data &photon,
                        random_dirty &random_obj)

{
    long i;
    int escape = 0;

#ifdef DEBUG_RT
    cout << "rt: io start; ";
    cout.flush();
#endif
    // initialize (and allocate if needed) output structure
    initialize_output(output, geometry);
#ifdef DEBUG_RT
    cout << "rt: io end; ";
    cout.flush();
#endif

    i = 0;
    int done_rt = 0;
    double total_stellar_weight_err = 0.0;
    double total_scattered_weight_err = 0.0;
    while ((i < geometry.n_photons) && (!done_rt))
    {
        if ((i > 0) && ((i % (geometry.n_photons / 20)) == 0) && (runinfo.rt_check_converged))
        {
            // compute total stellar weight uncertainty
            // first compute the uncertainty on the average weight of a single stellar
            // photon
            if (output.outputs[0].total_num_photons > 0)
                total_stellar_weight_err =
                    output.outputs[0].total_stellar_weight_x2 / output.outputs[0].total_num_photons -
                    pow(output.outputs[0].total_stellar_weight / output.outputs[0].total_num_photons, 2);
            if (total_stellar_weight_err > 0.0)
                total_stellar_weight_err = sqrt(total_stellar_weight_err / output.outputs[0].total_num_photons);
            else
                total_stellar_weight_err = 0.0;
            // this converts to a fractional uncertainty which is the same as the
            // fractional uncertainty on the total stellar weight
            total_stellar_weight_err /= output.outputs[0].total_stellar_weight / output.outputs[0].total_num_photons;

            // compute total scattered weight uncertainty
            // first compute the uncertainty on the average weight of a single
            // scattered photon
            if (output.outputs[0].total_num_photons > 0)
                total_scattered_weight_err =
                    output.outputs[0].total_scattered_weight_x2 / output.outputs[0].total_num_scattered_photons -
                    pow(output.outputs[0].total_scattered_weight / output.outputs[0].total_num_scattered_photons, 2);
            if (total_scattered_weight_err > 0.0)
                total_scattered_weight_err =
                    sqrt(total_scattered_weight_err / output.outputs[0].total_num_scattered_photons);
            else
                total_scattered_weight_err = 0.0;
            // this converts to a fractional uncertainty which is the same as the
            // fractional uncertainty on the total scattered weight
            total_scattered_weight_err /=
                output.outputs[0].total_scattered_weight / output.outputs[0].total_num_scattered_photons;

            if ((total_stellar_weight_err < runinfo.rt_converge_target) &&
                (total_scattered_weight_err < runinfo.rt_converge_target))
                done_rt = 1;
        }

        // print a status statement if asked
        if (runinfo.verbose >= 1)
        {
            if (((geometry.n_photons > 100) && ((i % (geometry.n_photons / 20)) == 0)) && (i > 0))
            {
                //|| ((i >= 25650) && (i <= 34200))) {
                cout << "current # = " << i;
                cout << " stel sl = " << output.outputs[0].total_stellar_weight / output.outputs[0].total_num_photons;
                cout << " scat sl = " << output.outputs[0].total_scattered_weight / output.outputs[0].total_num_photons;
                if (runinfo.rt_check_converged)
                    cout << "(" << total_scattered_weight_err << ")";
                cout << " ave 1st tau = " << output.outputs[0].ave_first_tau / output.outputs[0].total_num_photons;
                cout << endl;
            }
        }
        // start the photon at the star position
#ifdef DEBUG_RT
        cout << "rt: np begin; ";
        cout.flush();
#endif
#ifdef DEBUG_OUTRANGE
        try
        {
#endif
            photon.number = i;
            new_photon(photon, geometry, runinfo, random_obj);
#ifdef DEBUG_OUTRANGE
        }
        catch (std::out_of_range)
        {
            cout << "photon # = " << photon.number << endl;
            cout << "new_photon - out of range." << endl;
            exit(8);
        }
#endif
#ifdef DEBUG_RT
        cout << "rt: np end; ";
        cout.flush();
#endif
#ifdef DEBUG_RT
        if (photon.number > OUTNUM)
            cout << photon.number << " ";
        if (photon.number > OUTNUM)
            cout << "np1* done; ";
        cout.flush();
#endif

#ifdef SAVE_TRAJ
        cout << photon.number << " ";
        cout << photon.num_scat << " ";
        int k;
        for (k = 0; k < 3; k++)
            cout << photon.position[k] << " ";
        cout << " save_traj begin [pnum, nscat, pos(1-3)]";
        cout << endl;
#endif

        // find the first scattering site (forced)
#ifdef DEBUG_OUTRANGE
        try
        {
#endif

            escape = forced_first_scatter(geometry, photon, random_obj);

#ifdef DEBUG_OUTRANGE
        }
        catch (std::out_of_range)
        {
            cout << "photon # = " << photon.number << endl;
            cout << "force_first_scatter - out of range." << endl;
            exit(8);
        }
#endif
#ifdef DEBUG_RT
        if (photon.number >= OUTNUM)
            cout << "esc=" << escape << "; ";
        cout.flush();
        if (photon.number >= OUTNUM)
            cout << "ffs done; ";
        cout.flush();
#endif

#ifdef SAVE_TRAJ
        cout << photon.number << " ";
        cout << photon.num_scat << " ";
        for (k = 0; k < 3; k++)
            cout << photon.position[k] << " ";
        cout << " save_traj [pnum, nscat, pos(1-3)]";
        cout << endl;
#endif

        // classify stellar photon(s)
#ifdef DEBUG_RT
        if (photon.number >= OUTNUM)
            cout << "cstp start; ";
        cout.flush();
#endif
#ifdef DEBUG_OUTRANGE
        try
        {
#endif
            classify_stellar_photon(output, photon, geometry, runinfo, random_obj);
#ifdef DEBUG_OUTRANGE
        }
        catch (std::out_of_range)
        {
            cout << "photon # = " << photon.number << endl;
            cout << "classify_stellar_photon - out of range." << endl;
            exit(8);
        }
#endif
#ifdef DEBUG_RT
        if (photon.number >= OUTNUM)
            cout << "cstp done; ";
        cout.flush();
#endif

        // loop until photon escapes
        while (!escape)
        {
            // classify the scattered photon
#ifdef DEBUG_RT
            if (photon.number > OUTNUM)
                cout << "cscp begin; ";
            cout.flush();
#endif
#ifdef DEBUG_OUTRANGE
            try
            {
#endif
                classify_scattered_photon(output, photon, geometry, runinfo);
#ifdef DEBUG_OUTRANGE
            }
            catch (std::out_of_range)
            {
                cout << "photon # = " << photon.number << endl;
                cout << "num scat = " << photon.num_scat << endl;
                cout << "classify_scattered_photon - out of range." << endl;
                exit(8);
            }
#endif
#ifdef DEBUG_RT
            if (photon.number > OUTNUM)
                cout << "cscp done; ";
            cout.flush();
#endif

            // scatter the photon into a new direction
#ifdef DEBUG_OUTRANGE
            try
            {
#endif
                scatter_photon(geometry, photon, random_obj);
#ifdef DEBUG_OUTRANGE
            }
            catch (std::out_of_range)
            {
                cout << "photon # = " << photon.number << endl;
                cout << "num scat = " << photon.num_scat << endl;
                cout << "scatter_photon - out of range." << endl;
                exit(8);
            }
#endif
#ifdef DEBUG_RT
            if (photon.number > OUTNUM)
                cout << "sp done; ";
            cout.flush();
#endif

            // determine the next scattering site or
            // the position where the photon escapes the dust
#ifdef DEBUG_OUTRANGE
            try
            {
#endif
                escape = next_scatter(geometry, photon, random_obj);
#ifdef DEBUG_OUTRANGE
            }
            catch (std::out_of_range)
            {
                cout << "photon # = " << photon.number << endl;
                cout << "num scat = " << photon.num_scat << endl;
                cout << "next_scatter - out of range." << endl;
                exit(8);
            }
#endif
#ifdef DEBUG_RT
            if (photon.number > OUTNUM)
                cout << "ns done; ";
            cout.flush();
#endif

#ifdef SAVE_TRAJ
            cout << photon.number << " ";
            cout << photon.num_scat << " ";
            double rad = 0.0;
            for (k = 0; k < 3; k++)
            {
                cout << photon.position[k] << " ";
                rad += pow(photon.position[k], 2);
            }
            cout << escape << " ";
            cout << pow(rad, 0.5) << " ";
            cout << " save_traj [pnum, nscat, pos(1-3), escape, radius]";
            cout << endl;
#endif

#ifdef DEBUG_RT
            if (photon.number > OUTNUM)
                cout << endl;
#endif
        }

        i++;
    }
}

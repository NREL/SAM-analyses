#=
This code uses a branch of REopt.jl, which is deprecated.
To get the source code:
```sh
git clone https://github.com/NREL/REopt.git
git checkout run_sam_battery_stateful
```
Then add the REoptLite code to your Julia env
```julia
pkg> add your/path/to/REoptLite/  # that you just cloned
```

NOTE the SAM libraries in REoptLite for Mac do not support Apple chips.
=#
using JuMP, HiGHS, JSON, DelimitedFiles, CSV, DataFrames#, REopt


"""
Given errors in battery dispatch relative to forecast, update ElectricTariff variables for peak demand charge
which carry forward to the next MPC run
"""
function update_electric_tariff_demand!(results::Dict, inputs::Dict, start_index::Integer, next_start::Integer, batt_forecast_error::Vector{<:Real})
    monthly_peak = inputs["ElectricTariff"]["monthly_previous_peak_demands"]
    tou_peaks = inputs["ElectricTariff"]["tou_previous_peak_demands"]
    tou_timesteps = inputs["ElectricTariff"]["tou_demand_timesteps"]

    load_from_batt = results["ElectricStorage"]["to_load_series_kw"]
    load_from_grid = results["ElectricUtility"]["to_load_series_kw"] # increase if discharging shortfall

    n = size(batt_forecast_error)[1]

    battery_from_grid = zeros(n)
    if (haskey(results["ElectricUtility"], "to_battery_series_kw"))
        battery_from_grid = results["ElectricUtility"]["to_battery_series_kw"] # increase if PV shortfall for charging step 
    end
    j = 1
    # By the time we got here batt_ac has been adjusted for the load, but storage to load 
    while j <= n

        if batt_forecast_error[j] < 0 && load_from_batt[j] > 0
            extra_grid = min(-1.0 * batt_forecast_error[j], load_from_batt[j])
            load_from_grid[j] += extra_grid
        end
        if batt_forecast_error[j] > 0 && battery_from_grid[j] > 0
            battery_from_grid[j] -= batt_forecast_error[j]
            battery_from_grid[j] = max(battery_from_grid[j], 0)
        end

        j += 1
    end

    total_from_grid = load_from_grid .+ battery_from_grid

    start_month = REopt.month_of_hour(start_index)
    next_month = REopt.month_of_hour(next_start)

    if start_month == next_month
        for i in collect(1:size(tou_peaks)[1])
            for t in tou_timesteps[i]
                power = total_from_grid[t]
                if power > tou_peaks[i]
                    tou_peaks[i] = power
                end
            end
        end
        power = maximum(total_from_grid)
        if power > monthly_peak[1]
            monthly_peak[1] = power
        end
    else
        for i in collect(1:size(tou_peaks)[1])
            tou_peaks[i] = 0
        end
        monthly_peak[1] = 0
    end

end

"""
Sum relevant sections of results dictionary to get a single AC power flow for SAM use
"""
function get_ac_batt_power(results::Dict{String, Any})::Vector{Float64}
    pv_to_battery = results["PV"]["to_battery_series_kw"]
    battery_to_load = results["ElectricStorage"]["to_load_series_kw"]
    n = length(pv_to_battery)
    grid_to_battery = zeros(n)
    if (haskey(results["ElectricUtility"], "to_battery_series_kw"))
        grid_to_battery = results["ElectricUtility"]["to_battery_series_kw"]
    end
    return battery_to_load - pv_to_battery - grid_to_battery
end


# Need to adjust the dict objects here, since they'll be re-used post battery run to account for SOC errors
"""
Apply forecast errors to results dictionary:
 - Increase grid usage if PV is insufficent to cover it
 - Decrease battery discharging if load has decreased to prevent battery to grid
This function changes values in the results dictionary, these updated values will be used to update the utility rate
    accounting in later functions.
Positive pv_forecast_error and load_forecast_error mean that load was forecast to be greater than actual
indicating increased grid use.
"""
function adjust_ac_power_for_forecast!(results::Dict{String, Any}, ac_power::Vector{<:Real}, pv_forecast_error::Vector{<:Real}, load_forecast_error::Vector{<:Real})
    pv_to_battery = results["PV"]["to_battery_series_kw"]
    load_from_pv = results["PV"]["to_load_series_kw"]
    load_from_batt = results["ElectricStorage"]["to_load_series_kw"]
    load_from_grid = results["ElectricUtility"]["to_load_series_kw"] # increase if discharging shortfall

    # Sign convention matches grid use - positive indicates an increase in grid use
    net_error = pv_forecast_error - load_forecast_error

    batt_power_series = zeros(0)
    n = length(pv_to_battery)
    i = 1
    while i <= n
        batt_power = ac_power[i]
        pv_error = pv_forecast_error[i] # Assign the error to reduced battery charging or increased grid use
        load_error = load_forecast_error[i]
        if net_error[i] > 0 # increase in net load relative to forecast (w/o accounting for battery)
            
            # Discharging - nothing to do, PV and load functions below will take care of it                
            # Charging - reduce battery power if PV is less than forecast
            if batt_power < 0
                if pv_to_battery[i] > 0 && pv_forecast_error[i] > 0
                    charging_diff = max(0, pv_to_battery[i] - pv_error)
                    pv_to_battery[i] -= charging_diff
                    pv_error = max(0, pv_error - charging_diff) # Error reduces battery charging first
                    batt_power = max(0, ac_power[i] + charging_diff)
                end
            end

            
            if pv_error > 0 && load_from_pv[i] > 0
                extra_grid = min(pv_error, load_from_pv[i])
                pv_error = max(0, pv_error - extra_grid)  # pv error unaccounted for by reassigning load_from_pv
                load_from_pv[i] -= extra_grid 
                load_from_grid[i] += extra_grid
            end

            load_from_grid[i] += pv_error - load_error
            
        elseif net_error[i] < 0  # decrease in net load relative to forecast (w/o accounting for battery)
            # Discharging - check that there is sufficient load to absorb the battery power
            if batt_power > 0  # battery discharging
                # Prevent discharging to grid if insufficent load
                if load_from_batt[i] > 0 && load_error < 0  
                    discharging_diff = min(load_from_batt[i], -1.0* load_error)
                    load_error += discharging_diff  # Load error unaccounted for by reducing battery discharge
                    load_from_batt[i] -= discharging_diff
                    batt_power = min(0, batt_power - discharging_diff)
                end
            end
            #  Nothing to do for charging in this case, sufficient power is available

            load_from_grid[i] += pv_error - load_error
        end

        append!(batt_power_series, batt_power)
        i += 1
    end
end


function main(last_time_step = 8760)
    input_data = JSON.parsefile(joinpath(@__DIR__, "reopt_results", "reopt_results_outage_True_True_18.389_-66.0933_match_sam.json"))

    site_dict = input_data["inputs"]["Scenario"]["Site"]
    output_dict = input_data["outputs"]["Scenario"]["Site"]

    pv_outputs = output_dict["PV"]
    pv_capacity = pv_outputs["size_kw"]

    batt_dict = site_dict["Storage"]
    batt_outputs = output_dict["Storage"]
    batt_capacity = batt_outputs["size_kwh"]
    batt_power_kw = batt_outputs["size_kw"]
    grid_charge = batt_dict["canGridCharge"]
    inv_eff = batt_dict["inverter_efficiency_pct"]
    rec_eff = batt_dict["rectifier_efficiency_pct"]
    internal_eff = batt_dict["internal_efficiency_pct"]
    soc = batt_dict["soc_init_pct"]

    charge_eff = rec_eff * sqrt(internal_eff)
    discharge_eff = inv_eff * sqrt(internal_eff)

    tariff = REopt.URDBrate("5bfdc7925457a33744146c53", 2018)

    periods_per_month = floor(Int, size(tariff.tou_demand_rates)[1] / 12)

    # Hours
    start_index = 1
    horizon = 24
    interval = 1

    scenario_dict = Dict(
        "PV" => Dict("size_kw" => pv_capacity, "production_factor_series" => zeros(horizon)),
        "ElectricStorage" => Dict(
            "size_kw" => batt_power_kw, 
            "size_kwh" => batt_capacity, 
            "can_grid_charge" => grid_charge, 
            "soc_init_fraction" => soc, 
            "charge_efficiency" => charge_eff, 
            "discharge_efficiency" => discharge_eff
        ),
        "ElectricLoad" => Dict("loads_kw" => zeros(horizon)),
        "ElectricUtility" => Dict(),
        "ElectricTariff" => Dict(
            "net_metering" => false, 
            "export_rates" => zeros(horizon), 
            "monthly_previous_peak_demands" => [0.0], 
            "tou_previous_peak_demands" => zeros(periods_per_month)
        )
    )

    # TODO - get SAM battery SOC to match scenario_dict above
    # REopt.SAM_Battery class will resize appropriately given the scenarios above
    batt = REopt.SAM_Battery(joinpath(@__DIR__, "test_batt.json"), batt_capacity)

    actual_pv_df = CSV.read(joinpath(@__DIR__,"pv_production_actual.csv"), DataFrame) # 25 years of data
    forecast_pv_df = CSV.read(joinpath(@__DIR__, "pv_production_forecast.csv"), DataFrame) # 25 years of data

    actual_load_df = CSV.read(joinpath(@__DIR__, "weather_and_load", "san_juan_hospital_actual_load.csv"), DataFrame)
    forecast_load_df = CSV.read(joinpath(@__DIR__, "weather_and_load", "san_juan_hospital_forecast_load.csv"), DataFrame)

    outage_start = 4543
    outage_end = 4567

    forecast_output_powers = zeros(0)
    actual_output_powers = zeros(0)

    run_forecast =false

    while start_index < last_time_step + 1
        end_index = min(last_time_step, start_index + horizon - 1)
        rate_dict = REopt.get_subset_of_urdb(tariff, start_index, end_index)
        for (key, value) in rate_dict
            scenario_dict["ElectricTariff"][key] = value
        end
        # Forecast
    
        scenario_dict["ElectricStorage"]["soc_min_fraction"] = 0.7
        scenario_dict["ElectricTariff"]["export_rates"] = repeat([0.075], end_index - start_index + 1)

        """ SAM will deal with grid outage, none of the other algorithms get advance notice
        if start_index - outage_start > 0 && start_index <= outage_end
            scenario_dict["ElectricUtility"] = Dict()
            outage_ts = max(1, outage_start - start_index)
            outage_end_ts = max(1, outage_end - start_index)
            scenario_dict["ElectricUtility"]["outage_start_time_step"] = outage_ts
            scenario_dict["ElectricUtility"]["outage_end_time_step"] = outage_end_ts
        else
            delete!(scenario_dict, "ElectricUtility")
        end
        """

        if run_forecast
            forecast_pv = [pv < 0.0 ? 0.0 : pv for pv in forecast_pv_df[start_index:end_index, :Power]]
            forecast_load = forecast_load_df[start_index:end_index, :Load]

            actual_pv = [pv < 0.0 ? 0.0 : pv for pv in actual_pv_df[start_index:end_index, :Power]]
            actual_load = actual_load_df[start_index:end_index, :Load]

            forecast_pv[1] = actual_pv[1]
            forecast_load[1] = actual_load[1]

            scenario_dict["PV"]["production_factor_series"] = forecast_pv / scenario_dict["PV"]["size_kw"]
            scenario_dict["ElectricLoad"]["loads_kw"] = forecast_load
        else
            actual_pv = [pv < 0.0 ? 0.0 : pv for pv in actual_pv_df[start_index:end_index, :Power]]
            actual_load = actual_load_df[start_index:end_index, :Load]
            
            forecast_pv = actual_pv 
            forecast_load = actual_load

            scenario_dict["PV"]["production_factor_series"] = actual_pv / scenario_dict["PV"]["size_kw"]
            scenario_dict["ElectricLoad"]["loads_kw"] = actual_load
        end

        model = Model(HiGHS.Optimizer)
        results = REopt.run_mpc(model, scenario_dict)
        ac_batt_power = get_ac_batt_power(results)

        # In both variables below, positive means increased grid use
        pv_forecast_error = forecast_pv - actual_pv
        load_forecast_error = forecast_load - actual_load
        
        adjust_ac_power_for_forecast!(results, ac_batt_power, pv_forecast_error, load_forecast_error)

        reopt_batt_dc_power = REopt.get_batt_power_time_series(results, inv_eff, rec_eff)

        # Run the battery for number of steps in the interval to get the updated actual SOC
        actual_power = REopt.run_sam_battery(batt, reopt_batt_dc_power[1:interval])
        #print(actual_power)
        ac_actual_power = REopt.dc_to_ac_power(actual_power, inv_eff, rec_eff)
        REopt.update_mpc_from_batt_stateful(batt, scenario_dict)

        # Compare ac_actual power to ac
        # Positive - charged less than expected; negative - discharged less than expected
        batt_forecast_error = ac_batt_power[1:interval] - ac_actual_power
        
        update_electric_tariff_demand!(results, scenario_dict, start_index, start_index + interval, batt_forecast_error)
        
        # Want to save both of these timeseries for analysis
        if (interval < horizon)
            append!(forecast_output_powers, ac_batt_power[1:interval])
            append!(actual_output_powers, actual_power[1:interval])
        else
            append!(forecast_output_powers, ac_batt_power)
            append!(actual_output_powers, actual_power)
        end
        
        """
        print("Start time ", start_index, "\n")
        print("MPC results:\n", results)
        print("\n\n")
        print("Next inputs: \n", scenario_dict)
        print("\n\n")
        """
        start_index += interval
    end

    if run_forecast
        DelimitedFiles.writedlm(joinpath(@__DIR__,"output_powers_sj_hospital_day_ahead_24_1_forecast_ac_export_rate.csv"), forecast_output_powers, ',')
        DelimitedFiles.writedlm(joinpath(@__DIR__,"output_powers_sj_hospital_day_ahead_24_1_actual_dc_export_rate.csv"), actual_output_powers, ',')
    else
        DelimitedFiles.writedlm(joinpath(@__DIR__,"output_powers_sj_hospital_actual_24_1_ac_export_rate.csv"), forecast_output_powers, ',')
        DelimitedFiles.writedlm(joinpath(@__DIR__,"output_powers_sj_hospital_actual_24_1_dc_export_rate.csv"), actual_output_powers, ',')
    end
end

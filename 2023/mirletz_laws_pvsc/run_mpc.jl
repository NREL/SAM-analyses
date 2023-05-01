include(raw"C:\Users\bmirletz\source\repos\REoptLite\src\REoptLite.jl")
include(raw"C:\Users\bmirletz\source\repos\REoptLite\src\mpc\sam_battery.jl")

using Main.REoptLite, JuMP, Cbc, JSON, DelimitedFiles

function update_electric_tariff_demand(results::Dict, inputs::Dict, start_index::Integer, next_start::Integer)
    monthly_peak = inputs["ElectricTariff"]["monthly_previous_peak_demands"]
    tou_peaks = inputs["ElectricTariff"]["tou_previous_peak_demands"]
    tou_timesteps = inputs["ElectricTariff"]["tou_demand_timesteps"]

    load_from_grid = results["ElectricUtility"]["to_load_series_kw"]
    battery_from_grid = results["ElectricUtility"]["to_battery_series_kw"]
    total_from_grid = load_from_grid .+ battery_from_grid

    start_month = REoptLite.month_of_hour(start_index)
    next_month = REoptLite.month_of_hour(next_start)

    if start_month == next_month
        for i in collect(1:size(tou_peaks)[1])
            for t in tou_timesteps[i]
                power = total_from_grid[t]
                if power > tou_peaks[i]
                    tou_peaks[i] = power
                end
                if power > monthly_peak[1]
                    monthly_peak[1] = power
                end
            end
        end
    else
        for i in collect(1:size(tou_peaks)[1])
            tou_peaks[i] = 0
        end
        monthly_peak[1] = 0
    end

end

function get_ac_batt_power(results::Dict{String, Any})::Vector{Float64}
    pv_to_battery = results["PV"]["to_battery_series_kw"]
    grid_to_battery = results["ElectricUtility"]["to_battery_series_kw"]
    battery_to_load = results["Storage"]["to_load_series_kw"]

    batt_power_series = zeros(0)
    n = length(pv_to_battery)
    i = 1
    while i <= n
        charge = pv_to_battery[i] + grid_to_battery[i]
        batt_power = battery_to_load[i] - charge
        append!(batt_power_series, batt_power)
        i += 1
    end

    return batt_power_series
end

input_data = JSON.parsefile(raw"C:\Users\bmirletz\source\repos\mpc_analysis\reopt_results_32.74_-117.17_10_13_pv_only.json")

site_dict = input_data["inputs"]["Scenario"]["Site"]

pv_dict = site_dict["PV"]
pv_capacity = pv_dict["max_kw"] # Needs to change to process sized outputs
prod_factors = pv_dict["prod_factor_series_kw"]

batt_dict = site_dict["Storage"]
batt_capacity = batt_dict["max_kwh"] # Needs to change to process sized outputs
batt_power_kw = batt_dict["max_kw"]
grid_charge = batt_dict["canGridCharge"]
inv_eff = batt_dict["inverter_efficiency_pct"]
rec_eff = batt_dict["rectifier_efficiency_pct"]
internal_eff = batt_dict["internal_efficiency_pct"]
soc = batt_dict["soc_init_pct"]

charge_eff = rec_eff * sqrt(internal_eff)
discharge_eff = inv_eff * sqrt(internal_eff)

loads = site_dict["LoadProfile"]["loads_kw"]

generator = site_dict["Generator"]
run_during_outage = generator["generator_only_runs_during_grid_outage"]
gen_kw = generator["max_kw"]

tariff_label = site_dict["ElectricTariff"]["urdb_label"]
tariff = REoptLite.URDBrate("5cb743065457a321559b6ec4", 2018)

periods_per_month = floor(Int, size(tariff.tou_demand_rates)[1] / 12)

# Hours
global start_index = 1
horizon = 24
interval = 24

scenario_dict = Dict("PV" => Dict("size_kw" => pv_capacity, "prod_factor_series" => zeros(horizon)),
                    "Storage" => Dict("size_kw" => batt_power_kw, "size_kwh" => batt_capacity, "can_grid_charge" => grid_charge, "soc_init_pct" => soc, "charge_efficiency" => charge_eff, "discharge_efficiency" => discharge_eff),
                    "ElectricLoad" => Dict("loads_kw" => zeros(horizon)),
                    "ElectricTariff" => Dict("net_metering" => false, "export_rates" => zeros(horizon), "monthly_previous_peak_demands" => [0.0], "tou_previous_peak_demands" => zeros(periods_per_month)))

# TODO - get SAM battery SOC to match scenario_dict above
batt = SAM_Battery(raw"C:\Users\bmirletz\source\repos\mpc_analysis\test_batt.json", batt_capacity)

output_powers = zeros(0)

while start_index < 8761
    end_index = start_index + horizon - 1
    rate_dict = REoptLite.get_subset_of_urdb(tariff, start_index, end_index)
    for (key, value) in rate_dict
        scenario_dict["ElectricTariff"][key] = value
    end
    scenario_dict["PV"]["prod_factor_series"] = prod_factors[start_index:end_index]
    scenario_dict["ElectricLoad"]["loads_kw"] = loads[start_index:end_index]

    local model = Model(Cbc.Optimizer)
    results = REoptLite.run_mpc(model, scenario_dict)
    local ac_batt_power = get_ac_batt_power(results)
    local batt_power = get_batt_power_time_series(results, inv_eff, rec_eff)
    
    run_sam_battery(batt, batt_power)
    update_mpc_from_batt_stateful(batt, scenario_dict)
    
    update_electric_tariff_demand(results, scenario_dict, start_index, start_index + interval)
    
    append!(output_powers, ac_batt_power)
    
    print("Start time ", start_index, "\n")
    print("MPC results:\n", results)
    print("\n\n")
    print("Next inputs: \n", scenario_dict)
    print("\n\n")
    global start_index += interval
end

print(output_powers)
DelimitedFiles.writedlm(raw"C:\Users\bmirletz\source\repos\mpc_analysis\output_powers_sd_hospital_24_forecast.csv", output_powers, ',')
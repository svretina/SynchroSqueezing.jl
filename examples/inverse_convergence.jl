using SynchroSqueezing

fs = 1024.0
n = 1024
target_frequency = 64.0
target_bin = round(Int, target_frequency * n / fs) + 1

println("voices/octave  abs(R/Cψ - 1)  reliable range (Hz)")
for voices in (4, 8, 16, 32)
    plan = plan_cwt(n; fs, voices_per_octave = voices)
    error = abs(plan.synthesis_response[target_bin] / plan.admissibility - 1)
    println(rpad(voices, 15), rpad(error, 17), plan.reliable_frequency_range)
end

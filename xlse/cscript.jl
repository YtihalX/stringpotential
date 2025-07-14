using Libdl
using LinearAlgebra

# Load the shared library
const libscript = Libdl.dlopen(joinpath(@__DIR__, "build/linux/x86_64/release/libscript.so"))
using Plots
using LaTeXStrings

include("constants.jl")

function conshellT(E::Vector{Cdouble}, len, C::Vector{Cdouble}, pNgauss, Lambda, epsilon)
    @time otr = ccall(Libdl.dlsym(libscript, :onshellT), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Ptr{Cdouble}, Cuint, Cdouble, Cdouble), E, len, C, pNgauss, Lambda, epsilon)
    ot = transpose(copy(unsafe_wrap(Array, otr, (len, 4), own=false)))
    level = getEvec(C[1])
    vls = filter(e -> e > E[1] && e < E[end], level)
    vline(vls, s=:dash, c=:grey, label=L"$E_i$")
    vline!(delta, s=:dash, label="thresholds", lw=0.8)
    vline!([m_Xb11P], s=:dash, label=L"\chi_{b1}(1P)")
    vline!([m_Xb12P], s=:dash, label=L"\chi_{b1}(2P)")
    vline!([m_Xb13P], s=:dash, label=L"\chi_{b1}(3P)")
    plot!(E, abs.(ot[1, :]), label=L"$T_{11}$", dpi=300)
    plot!(E, abs.(ot[3,:]), label=L"$T_{21}$")
    plot!(E, abs.(ot[4, :]), label=L"$T_{22}$")
    plot!(E, abs.(ot[2, :]), label=L"$T_{12}$")
    # ylims!(0, upper)
    ylims!(0, 1e4)
    xlims!(E[1], E[end])
    savefig("onshellT.png")
    savefig("onshellT.pdf")
    cfree(reinterpret(Ptr{Cvoid}, otr))
    return ot
end

function conshellTV(E::Vector{Cdouble}, len, C::Vector{Cdouble}, pNgauss, Lambda, epsilon)
    @time otr = ccall(Libdl.dlsym(libscript, :onshellTV), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Ptr{Cdouble}, Cuint, Cdouble, Cdouble), E, len, C, pNgauss, Lambda, epsilon)
    ot = permutedims(copy(unsafe_wrap(Array, otr, (len, 4, 2), own=false)), (3, 2, 1))
    level = getEvec(C[1])
    vls = filter(e -> e > E[1] && e < E[end], level)
    vline(vls, s=:dash, c=:grey, label=L"$E_i$")
    vline!(delta, s=:dash, label="thresholds", lw=0.8)
    vline!([m_Xb11P], s=:dash, label=L"\chi_{b1}(1P)")
    vline!([m_Xb12P], s=:dash, label=L"\chi_{b1}(2P)")
    vline!([m_Xb13P], s=:dash, label=L"\chi_{b1}(3P)")
    plot!(E, abs.(ot[1, 1, :]) ./ 3 .- 700, label=L"$T_{11}$", dpi=300)
    # plot!(E, abs.(ot[1, 3,:]), label=L"$T_{21}$")
    # plot!(E, abs.(ot[1, 4, :]), label=L"$T_{22}$")
    # plot!(E, abs.(ot[1, 2, :]), label=L"$T_{12}$")
    plot!(E, abs.(ot[2, 1, :]), label=L"$V_{11}$", dpi=300)
    # plot!(E, abs.(ot[2, 3,:]), label=L"$V_{21}$")
    # plot!(E, abs.(ot[2, 4, :]), label=L"$V_{22}$")
    # plot!(E, abs.(ot[2, 2, :]), label=L"$V_{12}$")
    # ylims!(0, upper)
    ylims!(0, 500)
    xlims!(E[1], E[end])
    savefig("onshellTV.png")
    # savefig("onshellTV.pdf")
    cfree(reinterpret(Ptr{Cvoid}, otr))
    return ot
end

function conshellG(E::Vector{Cdouble}, len, C::Vector{Cdouble}, pNgauss, Lambda, epsilon)
    @time otr = ccall(Libdl.dlsym(libscript, :onshellG), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Ptr{Cdouble}, Cuint, Cdouble, Cdouble), E, len, C, pNgauss, Lambda, epsilon)
    ot = transpose(copy(unsafe_wrap(Array, otr, (len, 4), own=false)))
    upper = min(1e4, maximum(abs.(ot)))
    level = getEvec(C[1])
    vls = filter(e -> e > E[1] && e < E[end], level)
    vline(vls, s=:dash, c=:grey, label=L"$E_i$")
    vline!(delta, s=:dash, label="thresholds", lw=0.8)
    vline!([m_Xb11P], s=:dash, label=L"\chi_{b1}(1P)")
    vline!([m_Xb12P], s=:dash, label=L"\chi_{b1}(2P)")
    vline!([m_Xb13P], s=:dash, label=L"\chi_{b1}(3P)")
    plot!(E, abs.(ot[1, :]), label=L"$G_{11}$", dpi=300)
    # plot!(E, abs.(ot[3,:]), label=L"$G_{21}$")
    plot!(E, abs.(ot[4, :]), label=L"$G_{22}$")
    plot!(E, abs.(ot[2, :]), label=L"$G_{12}$")
    ylims!(0, upper)
    xlims!(E[1], E[end])
    # ylims!(0,10)
    savefig("onshellG.png")
    savefig("onshellG.pdf")
    cfree(reinterpret(Ptr{Cvoid}, otr))
    return ot
end

function conshellV(E::Vector{Cdouble}, len, C::Vector{Cdouble}, pNgauss, Lambda, epsilon)
    @time otr = ccall(Libdl.dlsym(libscript, :onshellV), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Ptr{Cdouble}, Cuint, Cdouble, Cdouble), E, len, C, pNgauss, Lambda, epsilon)
    ot = copy(transpose(unsafe_wrap(Array, otr, (len, 4), own=false)))
    level = getEvec(C[1])
    vls = filter(e -> e > E[1] && e < E[end], level)
    vline(vls, s=:dash, c=:grey, label=L"$E_i$")
    vline!(delta, s=:dash, label="thresholds", lw=0.8)
    vline!([m_Xb11P], s=:dash, label=L"\chi_{b1}(1P)")
    vline!([m_Xb12P], s=:dash, label=L"\chi_{b1}(2P)")
    vline!([m_Xb13P], s=:dash, label=L"\chi_{b1}(3P)")
    plot!(E, abs.(ot[1, :]), label=L"$V_{11}$", dpi=300)
    plot!(E, abs.(ot[3, :]), label=L"$V_{21}$")
    # plot!(E, abs.(ot[4, :]), label=L"$V_{22}$")
    plot!(E, abs.(ot[2, :]), label=L"$V_{12}$", dpi=300)
    xlims!(E[1], E[end])
    # ylims!(0, upper)
    ylims!(0, 5e3)
    savefig("onshellV.png")
    savefig("onshellV.pdf")
    cfree(reinterpret(Ptr{Cvoid}, otr))
    return ot
end

function detImVG(E::Vector{Cdouble}, len, C::Vector{Cdouble}, pNgauss, Lambda, epsilon)
    @time dtr = ccall(Libdl.dlsym(libscript, :Det), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Ptr{Cdouble}, Cuint, Cdouble, Cdouble), E, len, C, pNgauss, Lambda, epsilon)
    Det = copy(unsafe_wrap(Array, dtr, len, own=false))
    # plot(E, real.(Det), label=L"det($1-VG$)",dpi=300)
    level = getEvec(C[1])
    vls = filter(e -> e > E[1] && e < E[end], level)
    vline(vls, s=:dash, c=:grey, label=L"$E_i$")
    vline!(delta, s=:dash, label="thresholds", lw=0.8)
    vline!([m_Xb11P], s=:dash, label=L"\chi_{b1}(1P)")
    vline!([m_Xb12P], s=:dash, label=L"\chi_{b1}(2P)")
    vline!([m_Xb13P], s=:dash, label=L"\chi_{b1}(3P)")
    # vline!([m_pi + m_B - m_B_star], s=:dash, label=L"\pi")
    # vline!([m_pi], s=:dash, label=L"m_\pi")
    plot!(E, abs.(Det), label=L"|det($1-VG$)|", dpi=300)
    xlims!(E[1], E[end])
    ylims!(0,1e5)
    xlabel!("E/GeV")
    # ylims!(0, 1e9)
    savefig("det.png")
    savefig("det.pdf")
    cfree(reinterpret(Ptr{Cvoid}, dtr))
    return Det
end

function traceG(E::Vector{Cdouble}, len, C::Vector{Cdouble}, pNgauss, Lambda, epsilon)
    @time dtr = ccall(Libdl.dlsym(libscript, :traceG), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Ptr{Cdouble}, Cuint, Cdouble, Cdouble), E, len, C, pNgauss, Lambda, epsilon)
    trg = transpose(copy(unsafe_wrap(Array, dtr, (len, 2), own=false)))
    vline(delta, s=:dash, label="thresholds", lw=0.8)
    plot!(E, real.(trg[1, :]), label=L"real G_{11}", dpi=300)
    plot!(E, imag.(trg[1, :]), label=L"imag G_{11}", dpi=300)
    plot!(E, real.(trg[2, :]), label=L"real G_{22}", dpi=300)
    plot!(E, imag.(trg[2, :]), label=L"imag G_{22}", dpi=300)
    ylims!(-1.5, 0.5)
    savefig("trg.png")
    cfree(reinterpret(Ptr{Cvoid}, dtr))
    return trg
end

function lse_traceG(E::Vector{Cdouble}, C::Vector{Cdouble}, pNgauss, Lambda, epsilon)
    lse = ccall(Libdl.dlsym(libscript, :lse_malloc), Ptr{Cvoid}, (Csize_t, Cdouble, Cdouble), pNgauss, Lambda, epsilon)
    n = 2 * (pNgauss + 1)
    ch1 = []
    ch2 = []
    for e in E
        ccall(Libdl.dlsym(libscript, :lse_refresh), Cvoid, (Ptr{Cvoid}, ComplexF64, Ptr{Cdouble}, Cuint), lse, e, C, 3)
        ccall(Libdl.dlsym(libscript, :lse_gmat), Cint, (Ptr{Cvoid},), lse)
        ptr = ccall(Libdl.dlsym(libscript, :lse_get_g_data), Ptr{ComplexF64}, (Ptr{Cvoid},), lse)
        gg = copy(unsafe_wrap(Array, ptr, (n, n), own=false))
        push!(ch1, tr(gg[1:pNgauss+1, 1:pNgauss+1]))
        push!(ch2, tr(gg[pNgauss+2:end, pNgauss+2:end]))
    end
    plot(E, real.(ch1), dpi=300)
    plot!(E, imag.(ch1))
    plot!(E, real.(ch2))
    plot!(E, imag.(ch2))
    savefig("trg.png")
end

function Both(E::Vector{Cdouble}, len, C::Vector{Cdouble}, pNgauss, Lambda, epsilon)
    @time dtr = ccall(Libdl.dlsym(libscript, :Both), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Ptr{Cdouble}, Cuint, Cdouble, Cdouble), E, len, C, pNgauss, Lambda, epsilon)
    data = copy((unsafe_wrap(Array, dtr, (5, len), own=false)))
    ot = data[2:5, :]
    plot(E, abs.(ot[1, :]), label=L"$T_{11}$", dpi=300)
    plot!(E, abs.(ot[2, :]), label=L"$T_{12}$")
    # plot!(E, abs.(ot[3,:]), label=L"$T_{21}$")
    plot!(E, abs.(ot[4, :]), label=L"$T_{22}$")
    # ylims!(0,1e4)
    # savefig("onshellT.png")
    plot!(E, abs.(data[1, :]), label=L"|det($1-VG$)|", dpi=300)
    # ylims!(0,1.4)
    # ylims!(0, 1.2)
    level = getEvec(-1.01059)
    vls = filter(e -> e > E[1] && e < E[end], level)
    vline!(vls, s=:dash, c=:grey, label=L"$E_i$")
    ylims!(0, 1e4)
    savefig("det.png")
    savefig("det.pdf")
    cfree(reinterpret(Ptr{Cvoid}, dtr))
    return data
end

function Poles(Er::Vector{Cdouble}, rlen, Ei::Vector{Cdouble}, ilen, C::Vector{Cdouble}, pNgauss, Lambda, epsilon)
    @time dtr = ccall(Libdl.dlsym(libscript, :Poles), Ptr{ComplexF64}, (Ptr{Cdouble}, Csize_t, Ptr{Cdouble}, Csize_t, Ptr{Cdouble}, Csize_t, Cdouble, Cdouble), Er, rlen, Ei, ilen, C, pNgauss, Lambda, epsilon)
    data = copy(unsafe_wrap(Array, dtr, rlen * ilen * 4, own=false))
    data = transpose(reshape(data, (rlen * ilen, 4)))
    poles = [filter(!isnan, data[i, :]) for i in 1:4]
    cfree(reinterpret(Ptr{Cvoid}, dtr))
    # scatter(real.(poles[1]), imag.(poles[1]), dpi=300, label=L"RS_{--}", markersize=3, markerstrokewidth=0, legend=:outertopright)
    # scatter!(real.(poles[2]), imag.(poles[2]), dpi=300, label=L"RS_{+-}", markersize=3, markerstrokewidth=0)
    # scatter!(real.(poles[3]), imag.(poles[3]), dpi=300, label=L"RS_{-+}", markersize=3, markerstrokewidth=0)
    scatter(real.(poles[4]), imag.(poles[4]), dpi=300, label=L"RS_{++}", markersize=3, markerstrokewidth=0)
    level = getEvec(-1.01059)
    upper = -1e6
    lower = 1e6
    for i in 1:4
        upper = max(maximum(real.(poles[i])), upper)
        lower = min(minimum(real.(poles[i])), lower)
    end
    # println(level)
    vls = filter(e -> e > lower - 0.2 && e < upper + 0.2, level)
    vline!(vls, s=:dash, c=:grey, label=L"$E_i$")
    vline!(delta, s=:dash, label="thresholds", lw=0.8)
    vline!([m_Xb11P], s=:dash, label=L"\chi_{b1}(1P)")
    vline!([m_Xb12P], s=:dash, label=L"\chi_{b1}(2P)")
    vline!([m_Xb13P], s=:dash, label=L"\chi_{b1}(3P)")
    xlabel!("Re")
    ylabel!("Im")
    savefig("pole.png")
    savefig("pole.pdf")
    return poles
end

function minimize(C::Vector{Cdouble}, pNgauss, Lambda, epsilon)
    @time dtr = ccall(Libdl.dlsym(libscript, :Fit), Ptr{Cdouble}, (Ptr{Cdouble}, Csize_t, Cdouble, Cdouble), C, pNgauss, Lambda, epsilon)
    data = copy(unsafe_wrap(Array, dtr, 4, own=false))
    return data
end

function cost(c::Vector{Float64}, len, start::ComplexF64, pNgauss, Lambda, epsilon)
    @time dtr = ccall(Libdl.dlsym(libscript, :Cost), Ptr{Cdouble}, (Ptr{Cdouble}, Csize_t, ComplexF64, Csize_t, Cdouble, Cdouble), c, len, start, pNgauss, Lambda, epsilon)
    data = copy(unsafe_wrap(Array, dtr, len, own=false))
    cfree(reinterpret(Ptr{Cvoid}, dtr))
    plot(c, data, dpi=300, label="cost function")
    savefig("cost.png")
    savefig("cost.pdf")
end

function getEvec(V0)
    len = Ref{UInt}(0)
    ccall(Libdl.dlsym(libscript, :Ntower), Cvoid, (Ptr{Csize_t},), len)
    E = Array{Float64}(undef, len[])
    ccall(Libdl.dlsym(libscript, :Evec), Cvoid, (Ptr{Cdouble}, Cdouble), E, V0)
    return E
end

function cfree(ptr::Ptr{Cvoid})
    ccall(Libdl.dlsym(libscript, :Free), Cvoid, (Ptr{Cvoid},), ptr)
end

epsi = 1e-9
Lambda = 4
pNgauss = 40
data = Nothing
C = [-1.010589943548671, 0, -1.220749787118462, 0]
if "--poles" in ARGS
    Er = LinRange(-0.9, 0.0, 64)
    Ei = LinRange(-0.01, 0.01, 8)
    data = Poles(collect(Er), length(Er), collect(Ei), length(Ei), C, 64, Lambda, epsi)
    # cs = collect(LinRange(0.05, 1.4, 22))
    # po = []
    # for c in cs
    #     cxx = [c, 0, 0, 0]
    #     data = Poles(collect(Er), length(Er), collect(Ei), length(Ei), cxx, 64, Lambda, epsi)
    #     push!(po, minimum(real.(data[4])))
    # end
    # scatter(cs, po, dpi=300)
    # savefig("tmp.png")
end
onshellRange = LinRange(m_Xb11P - 0.3, delta[1] , 3000)

if "--onshellT" in ARGS
    # E = 1.48:0.00001:1.499
    E = onshellRange
    # C = [-4.015485e-01, -1.722080e+00, -1.854979e-01, -2.092185e+00]
    # C = [0.0, 0, 0, 0]
    # E = LinRange(-1.8, -1.6, 10000)
    # E = LinRange(-0.8001, -0.7999, 5000)
    # E = LinRange(0.5, 2, 5000)
    # E = 0.399999999:0.00000000001:0.400000001
    # E = 0:0.01:1.8
    # println("in julia")
    # println(collect(E))
    ot = conshellT(collect(E), length(E), C, pNgauss, Lambda, epsi)
    # og = conshellG(collect(E), length(E), 100, 4, epsi)
    # Det = detImVG(collect(E), length(E), 64, 4, epsi)
    # data = Both(collect(E), length(E), 64, 4, epsi)
    # E = getEvec()
end
if "--onshellG" in ARGS
    # E = 1.48:0.00001:1.499
    E = onshellRange
    # E = LinRange(-1.8, -1.6, 10000)
    # E = LinRange(-0.8001, -0.7999, 5000)
    # E = LinRange(0.5, 2, 5000)
    # E = 0.399999999:0.00000000001:0.400000001
    # E = 0:0.01:1.8
    # println("in julia")
    # println(collect(E))
    # ot = conshellG(collect(E), length(E), 100, 4, epsi)
    og = conshellG(collect(E), length(E), C, pNgauss, Lambda, epsi)
    # Det = detImVG(collect(E), length(E), 64, 4, epsi)
    # data = Both(collect(E), length(E), 64, 4, epsi)
    # E = getEvec()
end
if "--onshellV" in ARGS
    # E = 1.48:0.00001:1.499
    # C = [-4.015485e-01, -1.722080e+00, -1.854979e-01, -2.092185e+00]
    C = [0, 0.0, 0, 0]
    E = onshellRange
    # E = LinRange(-1.8, -1.6, 10000)
    # E = LinRange(-0.8001, -0.7999, 5000)
    # E = LinRange(0.5, 2, 5000)
    # E = 0.399999999:0.00000000001:0.400000001
    # E = 0:0.01:1.8
    ov = conshellV(collect(E), length(E), C, pNgauss, Lambda, epsi)
end
if "--onshellTV" in ARGS
    # E = 1.48:0.00001:1.499
    C = [-4.015485e-01, -1.722080e+00, -1.854979e-01, -2.092185e+00]
    E = onshellRange
    # E = LinRange(-1.8, -1.6, 10000)
    # E = LinRange(-0.8001, -0.7999, 5000)
    # E = LinRange(0.5, 2, 5000)
    # E = 0.399999999:0.00000000001:0.400000001
    # E = 0:0.01:1.8
    ov = conshellTV(collect(E), length(E), C, pNgauss, Lambda, epsi)
end

if "--traceG" in ARGS
    E = LinRange(-2, 1, 1000)
    data = traceG(collect(E), length(E), C, pNgauss, Lambda, epsi)
end

if "--Det" in ARGS
    # C = [-4.015485e-01, -1.722080e+00, -1.854979e-01, -2.092185e+00]
    # C = zeros(Float64, 4)
    E = onshellRange
    det = detImVG(collect(E), length(E), C, pNgauss, Lambda, epsi)
end

if "--Fit" in ARGS
    lse = ccall(dlsym(libscript, :lse_malloc), Ptr{Cvoid}, (Csize_t, Cdouble, Cdouble), pNgauss, Lambda, epsi)
    Er = LinRange(-0.98, 0.5, 64)
    Ei = LinRange(-0.01, 0.01, 16)
    # data = Poles(collect(Er), length(Er), collect(Ei), length(Ei), C, 64, Lambda, epsi)
    function deviation(c::Array{Float64})
        # c = [-1.01; c]
        tmp = ccall(dlsym(libscript, :lse_detImVG), ComplexF64, (Ptr{Cvoid}, ComplexF64, Ptr{Cdouble}, UInt), lse, m_Xb11P, c, 3)
        res = abs(tmp)
        tmp = ccall(dlsym(libscript, :lse_detImVG), ComplexF64, (Ptr{Cvoid}, ComplexF64, Ptr{Cdouble}, UInt), lse, m_Xb12P, c, 3)
        res += abs(tmp)
        tmp = ccall(dlsym(libscript, :lse_detImVG), ComplexF64, (Ptr{Cvoid}, ComplexF64, Ptr{Cdouble}, UInt), lse, m_Xb13P, c, 3)
        res += abs(tmp)
        return res
    end
    # deviation([0., 0, 0])
    using IMinuit
    fit = Minuit(deviation, (-1, 0.0, -1.3, 0))
    fit.strategy = 2
    # Fit(C, pNgauss, Lambda, epsi)
end

if "--cost" in ARGS
    c = LinRange(-2, 2, 5000)
    cost(collect(c), length(c), -0.7150819414927397 + 0im, pNgauss, Lambda, epsi)
end

if "--minimize" in ARGS
    # C[1] = -2
    C = [-2, 0.0, 0, 0]
    lse = ccall(dlsym(libscript, :lse_malloc), Ptr{Cvoid}, (Csize_t, Cdouble, Cdouble), pNgauss, Lambda, epsi)
    c = minimize(C, pNgauss, Lambda, epsi)
    println(ccall(dlsym(libscript, :lse_cost), Cdouble, (Ptr{Cvoid}, Ptr{Cdouble}, UInt), lse, c, 3))
    conshellT(collect(onshellRange), length(onshellRange), c, pNgauss, Lambda, epsi)
end

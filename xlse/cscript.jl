using Libdl
using Serialization
using Polynomials
using LinearAlgebra

# Load the shared library
const libscript = Libdl.dlopen(joinpath(@__DIR__, "build/linux/x86_64/release/libscript.so"))
using Plots
using LaTeXStrings
include("constants.jl")

xsqrtright(x) = imag(x) >= 0 ? sqrt(x + 0im) : -sqrt(x - 0im)
xsqrtup(x) = imag(x) >= 0 && real(x) < 0 ? -sqrt(x + 0im) : sqrt(x - 0im)
xsqrtleft(x) = sqrt(Complex(x + 0im))
# xsqrt(x) = imag(x) < 0 && real(x) < 0 ? -sqrt(x - 0im) : sqrt(x + 0im);

epsi = 1e-9
Lambda = 4.0
pNgauss = 64
data = Nothing
C = [-1.010589943548671, 0, -1.220749787118462, 0]
Erange = LinRange(-0.37, 0.5, 1000)
# Erange = LinRange(-0.1, 2, 1000)
# onshellRange = LinRange(-0.7, 0.6, 1000)
# onshellRange = LinRange(0., 0.190229863, 8000)

V(E, p, pprime) = ccall(dlsym(libscript, :V), ComplexF64, (Cdouble, ComplexF64, ComplexF64), E, p, pprime)
Vquad(E, p, pprime) = ccall(dlsym(libscript, :Vquad), ComplexF64, (Cdouble, ComplexF64, ComplexF64), E, p, pprime)
Gamma(E, p) = ccall(dlsym(libscript, :Gamma), ComplexF64, (ComplexF64, ComplexF64), E, p)
function integrand(e, p1, p2, x)
    A = p1^2 + p2^2 + m_pi^2
    B = 2 * p1 * p2
    C = (p1^2 + p2^2) / 2 / m_B - e - (m_B_star - m_B) - 1e-7im
    # D = 2 * m_B_star + (p1^2 + p2^2) / 2 / m_B_star - e - 1e-7im
    Eg = xsqrtleft(A - B * x)
    Da = 1 / (Eg + C)
    return Da / Eg
end

function quadrate(e, p1, p2)
    A = p1^2 + p2^2 + m_pi^2
    B = 2 * p1 * p2
    a = 1.3abs(A / B)
    x, w = gauss(64, 0, 1)
    xpath = [-1 .+ x .* (a * im); (-1 + a * im) .+ x .* 2; (1 + a * im) .+ x .* (-a * im)]
    wpath = [w .* (a * im); w .* 2; w .* (-a * im)]
    res = 0
    for i in eachindex(xpath)
        res += integrand(e, p1, p2, xpath[i]) * wpath[i]
    end
    return res
end

function EG(p1, p2, m0)
    A = p1^2 + p2^2 + m0^2
    B = 2 * p1 * p2
    x = LinRange(-1, 1, 100)
    return A .- B .* x
end
function onshell(E::Vector{Cdouble}, len, pNgauss, Lambda, epsilon)
    @time otr = ccall(Libdl.dlsym(libscript, :onshell), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Cuint, Cdouble, Cdouble), E, len, pNgauss, Lambda, epsilon)
    ot = transpose(copy(unsafe_wrap(Array, otr, (len, 2), own=false)))
    cfree(reinterpret(Ptr{Cvoid}, otr))
    plot(dpi=400)
    plot!(E, imag.(ot[1, :]))
    savefig("onshell.png")
    return ot
end

function imT(E::Vector{Cdouble}, len, C::Vector{Cdouble}, pNgauss, Lambda, epsilon)
    @time otr = ccall(Libdl.dlsym(libscript, :onshellT), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Ptr{Cdouble}, Cuint, Cdouble, Cdouble), E, len, C, pNgauss, Lambda, epsilon)
    ot = transpose(copy(unsafe_wrap(Array, otr, (len, 4), own=false)))
    invT = Array{ComplexF64}(undef, 4, len)
    function ρ(e, i)
        if e <= delta[i]
            return 0
        end
        k = sqrt(2 * mu[i] * (e - delta[i]))
        return mu[i] * k / 2 / pi
        e += m11 + m12
        s = e^2
        k = sqrt((s - (m[i, 1] + m[i, 2])^2) * (s - (m[i, 1] - m[i, 2])^2)) / 2 / e
        # return k
        return k / 8 / pi / e
    end
    for i in 1:len
        m = [ot[1, i] ot[2, i]; ot[3, i] ot[4, i]]
        if abs(det(m)) < 1e-8
            invT[1, i] = 0
            invT[2, i] = 0
            invT[3, i] = 0
            invT[4, i] = 0
            continue
        end
        invm = inv(m)
        invT[1, i] = invm[1, 1]
        invT[2, i] = invm[1, 2]
        invT[3, i] = invm[2, 1]
        invT[4, i] = invm[2, 2]
    end
    plot(dpi=400, legend=:topleft)
    plot!(E, imag.(invT[1, :]), label="Im " * L"T_{11}", lw=1.0, alpha=0.7)
    # plot!(E, imag.(invT[2, :]), label="Im "*L"T^{-1}_{12}")
    # plot!(E, imag.(invT[3, :]), label="Im "*L"T^{-1}_{21}")
    plot!(E, imag.(invT[4, :]), label="Im " * L"T_{22}", lw=1.0, alpha=0.7)
    plot!(E, ρ.(E, 1), label=L"\rho_1(E)\Theta(E-m_B-m_{B^*})", s=:dash, lw=1.5)
    plot!(E, ρ.(E, 2), label=L"\rho_2(E)\Theta(E-m_{B_s} - m_{B_s^*})", s=:dash, lw=1.5)
    vline!(delta, label="thresholds", s=:dash, c=:grey)
    ylims!(-0.6, 1.5)
    xlabel!("E/GeV")
    savefig("imT.png")
    savefig("iminvT.pdf")
    cfree(reinterpret(Ptr{Cvoid}, otr))
end

function imTsing(E::Vector{Cdouble}, len, C::Vector{Cdouble}, pNgauss, Lambda, epsilon)
    @time otr = ccall(Libdl.dlsym(libscript, :onshellT_single), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Ptr{Cdouble}, Cuint, Cdouble, Cdouble), E, len, C, pNgauss, Lambda, epsilon)
    ot = transpose(copy(unsafe_wrap(Array, otr, (len, 2), own=false)))
    invT = inv.(ot)
    function ρ(e, i)
        if e <= delta[i]
            return 0
        end
        k = sqrt(2 * mu[i] * (e - delta[i]))
        return mu[i] * k / 2 / pi
        e += m11 + m12
        s = e^2
        k = sqrt((s - (m[i, 1] + m[i, 2])^2) * (s - (m[i, 1] - m[i, 2])^2)) / 2 / e
        # return k
        return k / 8 / pi / e
    end
    p1 = xsqrtleft.(2 * mu[1] .* (E))
    p2 = 0.0006
    v1 = V.(E .+ (m11 + m12), p1, p2)
    plot(dpi=400, legend=:bottomright)
    plot!(E, imag.(invT[1, :]), label="Im " * L"T^{-1}_{11}", alpha=0.5, lw=1)
    # plot!(E, imag.(invT[2, :]), label="Im "*L"T^{-1}_{12}")
    # plot!(E, imag.(invT[3, :]), label="Im "*L"T^{-1}_{21}")
    # plot!(E, imag.(invT[4, :]), label="Im "*L"T^{-1}_{22}", alpha=0.5, lw=1, s=:dot)
    plot!(E, ρ.(E, 1), label=L"\rho_1(E)\Theta(E-m_B-m_{B^*})", alpha=0.8, s=:dash)
    plot!(E, ρ.(E, 2), label=L"\rho_2(E)\Theta(E-m_B-m_{B^*})", alpha=0.8, s=:dash)
    # plot!(E, real.(v1))
    # plot!(E, imag.(v1), label="imag")
    # plot!(E, imag.(invT[2, :]))
    # plot!(E, ρ.(E, 2), label=L"\rho_2(E)\Theta(E-m_{B_s} - m_{B_s^*})")
    vline!(delta, label="thresholds", s=:dash, c=:grey)
    # vline!([-m_eta^2/2/mu[1]], label=L"m_\eta")
    # ylims!(-2, 1.5)
    xlabel!("E/GeV")
    savefig("imTsing.png")
    savefig("iminvTsing.pdf")
    cfree(reinterpret(Ptr{Cvoid}, otr))
end

function conshellT(E::Vector{Cdouble}, len, C::Vector{Cdouble}, pNgauss, Lambda, epsilon)
    @time otr = ccall(Libdl.dlsym(libscript, :onshellT), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Ptr{Cdouble}, Cuint, Cdouble, Cdouble), E, len, C, pNgauss, Lambda, epsilon)
    yup = 1e2
    ot = permutedims(copy(unsafe_wrap(Array, otr, (len, 3, 3), own=false)), (3, 2, 1))
    plot(dpi=400, legend=:topleft)
    # level = getEvec(C[1])
    # vls = filter(e -> e > E[1] && e < E[end], level)
    # vline(vls, s=:dash, c=:grey, label=L"$E_i$")
    vline!(delta, s=:dash, label="thresholds", lw=0.8)
    vline!([m_Xb11P], s=:dash, label=L"\chi_{b1}(1P)")
    vline!([m_Xb12P], s=:dash, label=L"\chi_{b1}(2P)")
    vline!([m_Xb13P], s=:dash, label=L"\chi_{b1}(3P)")
    vline!([m_Xb14P], s=:dash, label=L"\chi_{b1}(4P)")
    plot!(E, abs.(ot[1, 1, :]), label=L"|$T_{11}$|", dpi=400)
    # plot!(E, abs.(ot[2, 1, :]), label=L"|$T_{21}$|")
    plot!(E, abs.(ot[2, 2, :]), label=L"|$T_{22}$|")
    plot!(E, abs.(ot[3, 3, :]), label=L"|$T_{33}$|")
    annotate!(0.01, -0.075 * yup, text(L"BB^*", 8))
    annotate!(delta[2] + 0.01, -0.075 * yup, text(L"B_sB_s^*", 8))
    # plot!(E, abs.(ot[2, :]), label=L"$T_{12}$")
    # ylims!(0, upper)
    xlabel!("E/GeV")
    ylabel!(L"|T_{\alpha\beta}|" * "/GeV")
    ylims!(0, yup)
    println(E[end])
    xlims!(E[1], E[end])
    savefig("onshellT.png")
    savefig("onshellT.pdf")
    cfree(reinterpret(Ptr{Cvoid}, otr))
    return ot
end

function conshellT_single(E::Vector{Cdouble}, len, C::Vector{Cdouble}, pNgauss, Lambda, epsilon)
    @time otr = ccall(Libdl.dlsym(libscript, :onshellT_single), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Ptr{Cdouble}, Cuint, Cdouble, Cdouble), E, len, C, pNgauss, Lambda, epsilon)
    ot = copy(unsafe_wrap(Array, otr, len, own=false))
    level = getEvec(C[1])
    vls = filter(e -> e > E[1] && e < E[end], level)
    # vline(vls, s=:dash, c=:grey, label=L"$E_i$")
    vline!(delta, s=:dash, label="thresholds", lw=0.8)
    # vline!([m_Xb11P], s=:dash, label=L"\chi_{b1}(1P)")
    # vline!([m_Xb12P], s=:dash, label=L"\chi_{b1}(2P)")
    # vline!([m_Xb13P], s=:dash, label=L"\chi_{b1}(3P)")
    plot!(E, abs.(ot), label=L"$T_{11}$", dpi=400)
    # vline!([m_pi-m_B_star], s=:dash, c=:grey, label=L"BB\pi")
    # plot!(E, abs.(ot[3, :]), label=L"$T_{21}$")
    # plot!(E, abs.(ot[4, :]), label=L"$T_{22}$")
    # plot!(E, abs.(ot[2, :]), label=L"$T_{12}$")
    # ylims!(0, upper)
    # ylims!(0, 3e3)
    println(E[end])
    xlims!(E[1], E[end])
    # ylims!(0, 50)
    xlabel!("E/GeV")
    savefig("onshellT_single.png")
    savefig("onshellT_single.pdf")
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
    plot!(E, abs.(ot[1, 1, :]) ./ 3 .- 700, label=L"$T_{11}$", dpi=400)
    # plot!(E, abs.(ot[1, 3,:]), label=L"$T_{21}$")
    # plot!(E, abs.(ot[1, 4, :]), label=L"$T_{22}$")
    # plot!(E, abs.(ot[1, 2, :]), label=L"$T_{12}$")
    plot!(E, abs.(ot[2, 1, :]), label=L"$V_{11}$", dpi=400)
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
    plot!(E, abs.(ot[1, :]), label=L"$G_{11}$", dpi=400)
    # plot!(E, imag.(ot[3,:]), label=L"$G_{21}$")
    plot!(E, abs.(ot[4, :]), label=L"$G_{22}$")
    plot!(E, abs.(ot[2, :]), label=L"$G_{12}$")
    # ylims!(0, upper)
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
    p = plot(dpi=400, layout=(1, 2))
    # vline(vls, s=:dash, c=:grey, label=L"$E_i$", dpi=400)
    vline!(p[1], delta, s=:dash, label="thresholds", lw=0.8)
    vline!(p[1], [m_Xb11P], s=:dash, label=L"\chi_{b1}(1P)")
    vline!(p[1], [m_Xb12P], s=:dash, label=L"\chi_{b1}(2P)")
    vline!(p[1], [m_Xb13P], s=:dash, label=L"\chi_{b1}(3P)")
    plot!(p[1], E, real.(ot[1, :]), label=L"real($V_{11}$)", dpi=400)
    xlabel!(p[1], "real, E/GeV")
    # plot!(E, abs.(ot[3, :]), label=L"$V_{21}$")
    # plot!(E, abs.(ot[4, :]), label=L"$V_{22}$")
    # plot!(E, abs.(ot[2, :]), label=L"$V_{12}$")
    xlims!(E[1], E[end])
    vline!(p[2], delta, s=:dash, label="thresholds", lw=0.8)
    vline!(p[2], [m_Xb11P], s=:dash, label=L"\chi_{b1}(1P)")
    vline!(p[2], [m_Xb12P], s=:dash, label=L"\chi_{b1}(2P)")
    vline!(p[2], [m_Xb13P], s=:dash, label=L"\chi_{b1}(3P)")
    plot!(p[2], E, imag.(ot[1, :]), label=L"imag($V_{11}$)", dpi=400)
    xlabel!(p[2], "imag, E/GeV")
    # plot!(E, abs.(ot[3, :]), label=L"$V_{21}$")
    # plot!(E, abs.(ot[4, :]), label=L"$V_{22}$")
    # plot!(E, abs.(ot[2, :]), label=L"$V_{12}$", dpi=400)
    # xlims!(E[1], E[end])
    # ylims!(0, 50)
    # ylims!(0, 5e3)
    savefig("onshellV.png")
    savefig("onshellV.pdf")
    cfree(reinterpret(Ptr{Cvoid}, otr))
    return ot
end

function detImVG(E::Vector{Cdouble}, len, C::Vector{Cdouble}, pNgauss, Lambda, epsilon)
    @time dtr = ccall(Libdl.dlsym(libscript, :Det), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Ptr{Cdouble}, Cuint, Cdouble, Cdouble), E, len, C, pNgauss, Lambda, epsilon)
    yup = 5e3
    Det = copy(unsafe_wrap(Array, dtr, len, own=false))
    # plot(E, real.(Det), label=L"det($1-VG$)",dpi=400)
    # level = getEvec(C[1])
    # vls = filter(e -> e > E[1] && e < E[end], level)
    # vline(vls, s=:dash, c=:grey, label=L"$E_i$")
    plot(dpi=400, legend=:topleft)
    vline!(delta, s=:dash, label="thresholds", lw=0.8)
    vline!([m_Xb11P], s=:dash, label=L"\chi_{b1}(1P)")
    vline!([m_Xb12P], s=:dash, label=L"\chi_{b1}(2P)")
    vline!([m_Xb13P], s=:dash, label=L"\chi_{b1}(3P)")
    vline!([m_Xb14P], s=:dash, label=L"\chi_{b1}(4P)")
    annotate!(0.01, -0.075 * yup, text(L"BB^*", 8))
    annotate!(delta[2] + 0.01, -0.075 * yup, text(L"B_sB_s^*", 8))
    annotate!(m_pi + m_B - m_B_star, -0.019 * yup, text(L"BB\pi", 8))
    vline!([m_pi + m_B - m_B_star], s=:dash, label=L"BB\pi")
    # vline!([m_pi + m_B - m_B_star], s=:dash, label=L"\pi")
    # vline!([m_pi], s=:dash, label=L"m_\pi")
    plot!(E, abs.(Det), label=L"|det($I+VG$)|", dpi=400)
    xlims!(E[1], E[end])
    ylims!(0, yup)
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
    plot!(E, real.(trg[1, :]), label=L"real G_{11}", dpi=400)
    plot!(E, imag.(trg[1, :]), label=L"imag G_{11}", dpi=400)
    plot!(E, real.(trg[2, :]), label=L"real G_{22}", dpi=400)
    plot!(E, imag.(trg[2, :]), label=L"imag G_{22}", dpi=400)
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
    plot(E, real.(ch1), dpi=400)
    plot!(E, imag.(ch1))
    plot!(E, real.(ch2))
    plot!(E, imag.(ch2))
    savefig("trg.png")
end

function Both(E::Vector{Cdouble}, len, C::Vector{Cdouble}, pNgauss, Lambda, epsilon)
    @time dtr = ccall(Libdl.dlsym(libscript, :Both), Ptr{ComplexF64}, (Ptr{Cdouble}, Cuint, Ptr{Cdouble}, Cuint, Cdouble, Cdouble), E, len, C, pNgauss, Lambda, epsilon)
    data = copy((unsafe_wrap(Array, dtr, (5, len), own=false)))
    ot = data[2:5, :]
    plot(E, abs.(ot[1, :]), label=L"$T_{11}$", dpi=400)
    # plot!(E, abs.(ot[2, :]), label=L"$T_{12}$")
    # plot!(E, abs.(ot[3,:]), label=L"$T_{21}$")
    # plot!(E, abs.(ot[4, :]), label=L"$T_{22}$")
    # ylims!(0,1e4)
    # savefig("onshellT.png")
    plot!(E, abs.(data[1, :]), label=L"|det($1-VG$)|", dpi=400)
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
    # scatter(real.(poles[1]), imag.(poles[1]), dpi=400, label=L"RS_{--}", markersize=3, markerstrokewidth=0, legend=:outertopright)
    # scatter!(real.(poles[2]), imag.(poles[2]), dpi=400, label=L"RS_{+-}", markersize=3, markerstrokewidth=0)
    # scatter!(real.(poles[3]), imag.(poles[3]), dpi=400, label=L"RS_{-+}", markersize=3, markerstrokewidth=0)
    scatter(real.(poles[4]), imag.(poles[4]), dpi=400, label=L"RS_{++}", markersize=3, markerstrokewidth=0)
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
    plot(c, data, dpi=400, label="cost function")
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
    # scatter(cs, po, dpi=400)
    # savefig("tmp.png")
end

if "--onshellT" in ARGS
    # E = 1.48:0.00001:1.499
    E = Erange
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
if "--onshellT_single" in ARGS
    # E = 1.48:0.00001:1.499
    E = Erange
    # C = [-4.015485e-01, -1.722080e+00, -1.854979e-01, -2.092185e+00]
    # C = [0.0, 0, 0, 0]
    # E = LinRange(-1.8, -1.6, 10000)
    # E = LinRange(-0.8001, -0.7999, 5000)
    # E = LinRange(0.5, 2, 5000)
    # E = 0.399999999:0.00000000001:0.400000001
    # E = 0:0.01:1.8
    # println("in julia")
    # println(collect(E))
    ot = conshellT_single(collect(E), length(E), C, pNgauss, Lambda, epsi)
    # og = conshellG(collect(E), length(E), 100, 4, epsi)
    # Det = detImVG(collect(E), length(E), 64, 4, epsi)
    # data = Both(collect(E), length(E), 64, 4, epsi)
    # E = getEvec()
end
if "--onshellG" in ARGS
    # E = 1.48:0.00001:1.499
    E = Erange
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
    E = Erange
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
    E = Erange
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
    E = Erange
    det = detImVG(collect(E), length(E), C, pNgauss, Lambda, epsi)
    # serialize("detwithoutrecoil.dat", det)
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
    conshellT(collect(Erange), length(Erange), c, pNgauss, Lambda, epsi)
end

z0(E, m, p1, p2, m0) = ((E - 2m - (p1^2 + p2^2) / (2m))^2 - m0^2 - p1^2 - p2^2) / (-2p2 * p1)
z0E(p1, p2, m0) = (p1^2 + p2^2 + m0^2) / (2p2 * p1)
if "--OME" in ARGS
    ome = ccall(dlsym(libscript, :ome_malloc), Ptr{Cvoid}, ())
    quad(E, p, pprime) = ccall(dlsym(libscript, :quad), ComplexF64, (ComplexF64, ComplexF64, ComplexF64), E, p, pprime)
    ana(E, p, pprime) = ccall(dlsym(libscript, :juliana), ComplexF64, (ComplexF64, ComplexF64, ComplexF64), E, p, pprime)
    Erange = LinRange(-0.801, 0.298, 4000)
    # Erange = LinRange(-0.801, 0.298, 1000)
    # Erange = LinRange(0.183, 0.184, 1000)
    prange = LinRange(0.00001, 2, 1000)
    p0 = 0.003525
    Esing = 0.1835 + m11 + m12
    # @time vsur = [abs(V(Esing, prange[i], prange[j])) for i in eachindex(prange), j in eachindex(prange)]
    # surface(prange, prange, vsur, dpi=400)
    # xlabel!("p")
    # ylabel!("p'")
    # zlabel!("V")
    # # @time vvec = [abs(V(Erange[i] + m11 + m12, prange[j], 0.003525)) for i in eachindex(Erange), j in eachindex(prange)]
    # # p = surface(Erange, prange, vvec, dpi=400)
    # # xlabel!(p, "E")
    # # ylabel!(p, "p")
    # zlims!(0, 80)
    # savefig("surface.png")
    pole1 = [z0(Erange[i] + m11 + m12, m_B, p0, p0, m_pi) for i in eachindex(Erange)]
    vvec = [abs(V(Erange[i] + m11 + m12, 0.003525, 0.003525)) for i in eachindex(Erange)]
    scatter(Erange, pole1, dpi=400, label="singularity position", markersize=1.5, markerstrokewidth=0)
    plot!(Erange, vvec, dpi=400, label="abs(numerical V)")
    hline!([-0.3, 0.3], label="box")
    xlabel!("E/GeV")
    # ylims!(-0.5, maximum(vvec) * 1.1)
    ylims!(-0.5, 5.5)
    savefig("pole1.png")
end

function nonana(E, p, m1, m2, m0)::Vector{ComplexF64}
    M = m1 + m2
    poly = Polynomial([E^2 + M^2 + p^4 / 4 / m2 / m2 - 2E * M - E * p^2 / m2 + M * p^2 / m2 - p^2 - m0^2, 2p, (M - E) / m1 + p^2 / 2 / m1 / m2 - 1, 0, 1 / 4 / m1 / m1,])
    rts = roots(poly)
    for i in eachindex(rts)
        if abs(imag(rts[i])) < 1e-6
            if E - M - real(rts[i])^2 / 2 / m1 - p^2 / 2 / m2 < 0
                rts[i] = NaN
            end
        end
    end
    return rts
end

function getV(E)
    dtr = ccall(dlsym(libscript, :getV), Ptr{ComplexF64}, (Cdouble, Csize_t, Cdouble, Cdouble), E, pNgauss, Lambda, epsi)
    V = copy(transpose(unsafe_wrap(Array, dtr, (2 * pNgauss + 2, 2 * pNgauss + 2), own=false)))
    cfree(reinterpret(Ptr{Cvoid}, dtr))
    return V
end

if "--analyticity" in ARGS
    prange = LinRange(0.01, 2, 10000)
    # prange = [0.9]
    # prange = vcat(LinRange(0.01 + 0.0001im, 0.9 + 0.01im, 600), LinRange(0.9 + 0.01im, 0.9 + 0.16im, 300), LinRange(0.9 + 0.16im, 2 + 0.16im, 600), LinRange(2 + 0.16im, 2, 300))
    # Erange = LinRange(-0.4, 0.0, 20)
    # Erange = [0.6]
    @time nonanalyticity = [nonana(E + m_B + m_B_star, p, m_B, m_B, m_pi) for E in Erange, p in prange]
    z = Array{ComplexF64}(undef, size(Erange)[1], 4 * size(prange)[1])
    for i in eachindex(Erange)
        z[i, :] = vcat(nonanalyticity[i, :]...)
    end
    # Create plot
    p = plot(legend=false, xlabel="Re", ylabel="Im", title="Domain of non-analyticity", dpi=400)

    # Plot each group with unique color
    for i in 1:size(z, 1)
        re = real.(z[i, :])
        im = imag.(z[i, :])
        scatter!(p, re, im,
            # markersize=0.5,
            markerstrokewidth=0,
            color=i,  # Unique color per group
            alpha=0.6) # Semi-transparent for overlap visibility
    end
    # scatter!(p, real.(prange), imag.(prange), markersize=0.5, markerstrokewidth=0, label="path", legend=true)


    xlims!(-1, 3)
    # ylims!(-0.05, 0.05)
    @time savefig(p, "complex.png")

end

if "--nroots" in ARGS
    prange = LinRange(0.01, 2, 1000)
    E = 0.8
    @time nonanalyticity = [nonana(E + m_B + m_B_star, p, m_B, m_B, m_pi) for p in prange]
    all_points = vcat(nonanalyticity...)
    for sol in nonanalyticity
        count = 0
        for r in sol
            if real(r) > 0 && real(r) <= 2 && abs(imag(r)) < 1e-6
                count += 1
            end
        end
        if count > 1
            print(count)
        end
    end
    # x_min, x_max = minimum(real.(all_points)), maximum(real.(all_points))
    # y_min, y_max = minimum(imag.(all_points)), maximum(imag.(all_points))
    #
    # # Add some padding to the limits for better visuals
    # padding_x = (x_max - x_min) * 0.1
    # padding_y = (y_max - y_min) * 0.1
    # xlims_with_padding = (x_min - padding_x, x_max + padding_x)
    # ylims_with_padding = (y_min - padding_y, y_max + padding_y)
    #
    #
    # # --- 3. Create the Animation ---
    # # We use the @gif macro to loop through each element of our data array.
    # println("Generating animation... this may take a moment.")
    #
    # data = nonanalyticity
    # @gif for i in eachindex(nonanalyticity)
    #
    #     current_points = data[i]
    #
    #     # To connect the four points and form a closed shape, we append the first point to the end
    #     shape_points = [current_points..., current_points[1]]
    #
    #     # Extract the real (x) and imaginary (y) coordinates for plotting
    #     x_coords = real.(shape_points)
    #     y_coords = imag.(shape_points)
    #
    #     # Create the plot for the current frame
    #     plot(x_coords, y_coords,
    #         xlabel="Real Part",
    #         ylabel="Imaginary Part",
    #         title="Frame $i of $(length(data))",
    #
    #         # Use the pre-calculated limits
    #         xlims=xlims_with_padding,
    #         ylims=ylims_with_padding,
    #
    #         # Ensure a square aspect ratio, crucial for complex planes
    #         aspect_ratio=:equal,
    #
    #         # Use seriestype=:shape to fill the polygon
    #         seriestype=:shape,
    #         fillalpha=0.3,
    #
    #         # Also show the individual points clearly
    #         marker=(:circle, 5, 0.8, :red, stroke(0)),
    #
    #         # Style the connecting line
    #         linewidth=2,
    #         linecolor=:black,
    #
    #         # Disable the legend as it's not needed here
    #         legend=false
    #     )
    #
    # end fps = 30
    #
    # println("Animation saved as 'anim.gif'")
end

if "--imT" in ARGS
    imT(collect(Erange), length(Erange), C, pNgauss, Lambda, epsi)

end

if "--imTsing" in ARGS
    imTsing(collect(Erange), length(Erange), C, pNgauss, Lambda, epsi)
end

o1(p1) = 2 * m_B + (p1^2 + p1^2) / (2m_B)
o2(p1) = 2 * m_B_star + (p1^2 + p1^2) / (2m_B_star)
if "--delt" in ARGS
    E = Erange
    k = sqrt.(Complex.(2 * mu[1] .* E))
    delt = Array{ComplexF64}(undef, length(E))
    for i in eachindex(E)
        delt[i] = V(E[i] + m11 + m12, k[i], k[i])
    end
    plot(dpi=400)
    # plot!(E, real.(delt), label="real")
    A = 2 .* k .* k .+ m_pi^2
    B = 2 .* k .* k
    C = o1.(k) .- (E .+ m11 .+ m12)
    D = o2.(k) .- (E .+ m11 .+ m12)
    a = sqrt.(A .- B)
    b = sqrt.(A .+ B)
    # plot!(E, imag.(log.(C .+ D)))
    plot!(E, imag.(delt), label="imag", s=:dash)
    savefig("delt.png")
end

if "--V3d" in ARGS
    using QuadGK
    # V(E, p, pprime)::ComplexF64 = ccall(dlsym(libscript, :V), ComplexF64, (Cdouble, ComplexF64, ComplexF64), E, p, pprime)
    pNgauss = 256
    x, w = gauss(pNgauss, 0, Lambda)
    @time dtr = ccall(Libdl.dlsym(libscript, :V3d), Ptr{ComplexF64}, (Cdouble, Cuint, Cdouble, Cdouble), m_pi, pNgauss, Lambda, epsi)
    n = 2 * pNgauss + 2
    data = transpose(copy(unsafe_wrap(Array, dtr, (n, n), own=false)))
    # data = [abs(V(m_pi + m11 + m12, x[i], x[j])) for i in eachindex(x), j in eachindex(x)]
    surface(x, x, abs.(data[1:pNgauss, 1:pNgauss]), dpi=400, camera=(45, 30))
    xlabel!("p")
    ylabel!("p'")
    zlabel!("abs(V)")
    zlims!(0, 3000)
    savefig("surface.png")
end

if "--quadrature" in ARGS
    using QuadGK
    xim, wim = gauss(20, 0, 1)
    xre, wre = gauss(40, 0, 1)
    x = [-1 .+ xim .* 0.2im; (-1 + 0.2im) .+ xre .* 2; (1 + 0.2im) .+ xim .* -0.2im]
    w = [wim .* 0.2im; wre .* 2; wim .* -0.2im]
    A = 0.053938
    B = 0.001367
    C = -0.234736 - 0.000001im
    f(x) = 1 / xsqrtright(A - B * x) / (xsqrtright(A - B * x) + C)
    dtr = ccall(Libdl.dlsym(libscript, :getIntegrand), Ptr{ComplexF64}, (),)
    data = copy(unsafe_wrap(Array, dtr, 80, own=false))
    cfree(reinterpret(Ptr{Cvoid}, dtr))
    t = -1
    println("t: $t")
    println("xsqrt: $(xsqrtright(A-B*t))")
    println("xsqrt+C: $(xsqrtright(A-B*t)+C)")
    t = -1 + 0.00000000001im
    println("t: $t")
    println("xsqrt: $(xsqrtright(A-B*t))")
    println("xsqrt+C: $(xsqrtright(A-B*t)+C)")
end

if "--cut" in ARGS
    using QuadGK
    E = LinRange(m11 + m12 - 0.3, m11 + m12 + 0.3, 500)
    p1 = xsqrtleft.(2 * mu[1] .* (E .- (m11 + m12)))
    p2 = 1
    # vana = Vquad.(E, p1, p2)
    v1 = V.(E, p1, p1)
    using QuadGK
    # varray = Array{ComplexF64}(undef, length(E))
    # for i in eachindex(varray)
    #     varray[i] = quadgk(x -> integrand(x, E[i]), -1, 1)[1]
    # end
    p = plot(layout=(1, 2), dpi=400, size=(1000, 500))
    # plot!(p[1], E, real.(vana), label="Re[quadrature]")
    # plot!(p[2], E, imag.(vana), label="Im[quadrature]")
    plot!(p[1], E, real.(v1), label="Re[analytical expression]")
    plot!(p[2], E, imag.(v1), label="Im[analytical expression]")
    vline!(p[1], [m11 + m12], label="threshold", c=:grey, s=:dash)
    vline!(p[2], [m11 + m12], label="threshold", c=:grey, s=:dash)
    # vline!(p[1], label=L"2\mu E + p^2 + m_{pi}^2 = 0", [(p2^2 + m_pi^2) / (-2mu[1]) + m11 + m12])
    # vline!(p[2], label=L"2\mu E + p^2 + m_{pi}^2 = 0", [(p2^2 + m_pi^2) / (-2mu[1]) + m11 + m12])
    savefig("cut.png")
    savefig("cut.pdf")
end


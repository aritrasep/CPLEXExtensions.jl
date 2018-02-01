###############################################################################
#                                                                             #
#  This file is part of the julia module for Multi Objective Optimization     #
#  (c) Copyright 2017 by Aritra Pal                                           #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     # 
# copy of this software and associated documentation files (the "Software"),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
# Every publication and presentation for which work based on the Program or   #
# its output has been used must contain an appropriate citation and           #
# acknowledgment of the author(s) of the Program.                             #
###############################################################################

module CPLEXExtensions

using Modof, CPLEX

#####################################################################
# Creating a CPLEX Model                                            #
#####################################################################

@inbounds function cplex_model{T<:MOOInstance}(instance::T, obj_ind::Int64, log_file_name::Union{String, Void}=nothing)
    env = CPLEX.Env()
    if typeof(instance) == MOBPInstance
        model = CPLEX.cplex_model(env, f = vec(instance.c[obj_ind, :]), lb = zeros(size(instance.c)[2]), ub = ones(size(instance.c)[2]))
    else
        model = CPLEX.cplex_model(env, f = vec(instance.c[obj_ind, :]), lb = instance.v_lb, ub = instance.v_ub)
    end
    if log_file_name != nothing
        CPLEX.set_logfile(model.env, log_file_name)
    end
    CPLEX.set_param!(model.env, "CPX_PARAM_THREADS", 1)
    if typeof(instance) == MOLPInstance
        #CPLEX.set_param!(model.env, 1016, instance.ϵ)
    else
        CPLEX.set_param!(model.env, "CPX_PARAM_EPGAP", 1.0e-6)
    end
    sense, rhs = Char[], Float64[]
    for i in 1:length(instance.cons_lb)
        if instance.cons_lb[i] != -Inf
            if instance.cons_ub[i] != Inf
                push!(sense, '=')
            else
                push!(sense, '>')
            end
            push!(rhs, instance.cons_lb[i])
        else
            push!(sense, '<')
            push!(rhs, instance.cons_ub[i])
        end
    end
    CPLEX.add_constrs!(model, instance.A, sense, rhs)
    if typeof(instance) == MOLPInstance
        return model
    end
    var_types = fill('I', size(instance.A)[2])
    if typeof(instance) == MOIPInstance || typeof(instance) == MOBPInstance
    else
        for i in 1:length(instance.var_types)
            if instance.var_types[i] == :Cont
                var_types[i] = 'C'
            end
        end
    end
    CPLEX.set_vartype!(model, var_types)
    model
end

@inbounds function cplex_model{T<:BOOInstance}(instance::T, obj_ind::Int64, log_file_name::Union{String, Void}=nothing)
    env = CPLEX.Env()
    c = obj_ind==1 ? vec(instance.c1) : vec(instance.c2)
    if typeof(instance) == BOBPInstance
        model = CPLEX.cplex_model(env, f = c, lb = zeros(length(c)), ub = ones(length(c)))
    else
        model = CPLEX.cplex_model(env, f = c, lb = instance.v_lb, ub = instance.v_ub)
    end
    if log_file_name != nothing
        CPLEX.set_logfile(model.env, log_file_name)
    end
    CPLEX.set_param!(model.env, "CPX_PARAM_THREADS", 1)
    if typeof(instance) == BOLPInstance
        #CPLEX.set_param!(model.env, 1016, instance.ϵ)
    else
        CPLEX.set_param!(model.env, "CPX_PARAM_EPGAP", 1e-6)
    end
    sense, rhs = Char[], Float64[]
    for i in 1:length(instance.cons_lb)
        if instance.cons_lb[i] != -Inf
            if instance.cons_ub[i] != Inf
                push!(sense, '=')
            else
                push!(sense, '>')
            end
            push!(rhs, instance.cons_lb[i])
        else
            push!(sense, '<')
            push!(rhs, instance.cons_ub[i])
        end
    end
    CPLEX.add_constrs!(model, instance.A, sense, rhs)
    if typeof(instance) == BOLPInstance
        return model
    end
    var_types = fill('I', size(instance.A)[2])
    if typeof(instance) == BOIPInstance || typeof(instance) == BOBPInstance
    else
        for i in 1:length(instance.var_types)
            if instance.var_types[i] == :Cont
                var_types[i] = 'C'
            end
        end
    end
    CPLEX.set_vartype!(model, var_types)
    model
end

#####################################################################
# Modifying the Constraint Matrix                                   #
#####################################################################

@inbounds del_constrs!{T<:Signed}(model::CPLEX.Model, ind::T) = del_constrs!(model, ind, ind)

@inbounds function del_constrs!{T<:Signed}(model::CPLEX.Model, start_ind::T, end_ind::T)
    stat = CPLEX.@cpx_ccall(delrows, Cint, (
                        Ptr{Void},
                        Ptr{Void},
                        Cint,
                        Cint
                        ),
                        model.env.ptr, model.lp, convert(Cint, start_ind-1), convert(Cint, end_ind-1))
    if stat != 0
        throw(CplexError(model.env, stat))
    end
end

@inbounds chg_coeffs!{T<:Real, S<:Real}(model::CPLEX.Model, rowlist::T, collist::T, vallist::S) = chg_coeffs!(model, Cint[rowlist], Cint[collist], vallist)
@inbounds chg_coeffs!{T<:Real, S<:Real}(model::CPLEX.Model, rowlist::Vector{T}, collist::Vector{T}, vallist::Vector{S}) = chg_coeffs!(model, convert(Vector{Cint}, rowlist), convert(Vector{Cint}, collist), CPLEX.fvec(vallist))

@inbounds function chg_coeffs!(model::CPLEX.Model, rowlist::Vector{Cint}, collist::Vector{Cint}, vallist::CPLEX.FVec)

    (length(rowlist) == length(collist) == length(vallist)) || error("Inconsistent argument dimensions.")
    stat = CPLEX.@cpx_ccall(chgcoeflist, Cint, (
                        Ptr{Void},
                        Ptr{Void},
                        Cint,
                     	Ptr{Cint},
                     	Ptr{Cint},
                     	Ptr{Float64}
                        ),
                        model.env.ptr, model.lp, convert(Cint, length(rowlist)), rowlist-Cint(1), collist-Cint(1), vallist)
    if stat != 0
        throw(CplexError(model.env, stat))
    end
end

@inbounds function set_rhs!{T<:Signed, S<:AbstractFloat}(model::CPLEX.Model, inds::Vector{T}, vals::Vector{S})
    stat = CPLEX.@cpx_ccall(chgrhs, Cint, (
                      Ptr{Void},
                      Ptr{Void},
                      Cint,
                      Ptr{Cint},
                      Ptr{Cdouble}
                      ),
                      model.env.ptr, model.lp, convert(Cint, length(inds)), inds - 1, vals)
    if stat != 0
        throw(CplexError(model.env, stat))
    end
end

#####################################################################
# Modifying an arbitary coefficient in the model                    #
#####################################################################

@inbounds chg_coeff_of_obj!{T<:Signed, S<:AbstractFloat}(model::CPLEX.Model, col::T, val::S) = chg_coeff!(model, 0, col, val)

@inbounds chg_coeff_of_rhs!{T<:Signed, S<:AbstractFloat}(model::CPLEX.Model, row::T, val::S) = chg_coeff!(model, convert(Cint, row), convert(Cint, 0), val)

@inbounds function chg_coeff!{T<:Signed, S<:AbstractFloat}(model::CPLEX.Model, row::T, col::T, val::S)

    stat = CPLEX.@cpx_ccall(chgcoef, Cint, (
                        Ptr{Void},
                        Ptr{Void},
                        Cint,
                     	Cint,
                     	Cdouble
                        ),
                        model.env.ptr, model.lp, row-1, col-1, val)
    if stat != 0
        throw(CplexError(model.env, stat))
    end
end

@inbounds function get_rhs_coef(model::CPLEX.Model, ind::Int64)
    rhs = Vector{Cdouble}(1)
    stat = CPLEX.@cpx_ccall(getrhs, Cint, (
                      Ptr{Void},
                      Ptr{Void},
                      Ptr{Cdouble},
                      Cint,
                      Cint
                      ),
                      model.env.ptr, model.lp, rhs, convert(Cint, ind)-1, convert(Cint, ind)-1)
    if stat != 0
        throw(CplexError(model.env, stat))
    end
    return rhs
end

#####################################################################
# MIP Warm Starts                                                   #
#####################################################################

@inbounds repair_an_infeasible_warm_start!{T<:Signed}(model::CPLEX.Model, x::BOPSolution, repair_tries::T=0) = repair_an_infeasible_warm_start!(model, x.vars, repair_tries)

@inbounds function repair_an_infeasible_warm_start!{T<:Signed}(model::CPLEX.Model, x::Vector{Float64}, repair_tries::T=0)
    CPLEX.set_param!(model.env, 2067, repair_tries)
    CPLEX.set_warm_start!(model, x, CPLEX.CPX_MIPSTART_REPAIR)
end

@inbounds add_a_feasible_warm_start!(model::CPLEX.Model, x::MOPSolution) = add_a_feasible_warm_start!(model, x.vars)
@inbounds add_a_feasible_warm_start!(model::CPLEX.Model, x::BOPSolution) = add_a_feasible_warm_start!(model, x.vars)

@inbounds add_a_feasible_warm_start!(model::CPLEX.Model, x::Vector{Float64}) = CPLEX.set_warm_start!(model, x, CPLEX.CPX_MIPSTART_NOCHECK)

@inbounds function del_all_warm_starts!(model::CPLEX.Model)
    try
        while true
            del_warm_start!(model, 1)
        end
    catch
    end
end

@inbounds del_warm_start!{T<:Signed}(model::CPLEX.Model, ind::T) = del_warm_starts!(model, convert(Cint, ind), convert(Cint, ind))

@inbounds del_warm_starts!{T<:Signed}(model::CPLEX.Model, start_ind::T, end_ind::T) = del_warm_starts!(model, convert(Cint, start_ind), convert(Cint, end_ind))

@inbounds function del_warm_starts!(model::CPLEX.Model, start_ind::Cint, end_ind::Cint)
    stat = CPLEX.@cpx_ccall(delmipstarts, Cint, (
                      Ptr{Void},
                      Ptr{Void},
                      Cint,
                      Cint
                      ),
                      model.env.ptr, model.lp, start_ind-1, end_ind-1)
    if stat != 0
        throw(CplexError(model.env, stat))
    end
end

#####################################################################
# Exporting Functions                                               #
#####################################################################

export cplex_model
export del_constrs!, chg_coeffs!, set_rhs!, chg_coeff_of_obj!, chg_coeff_of_rhs!, chg_coeff!, get_rhs_coef
export repair_an_infeasible_warm_start!, add_a_feasible_warm_start!, del_all_warm_starts!, del_warm_start!, del_warm_starts!

end

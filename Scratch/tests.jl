using Test
      
@test begin
        include("assignment_tests.jl")
        isapprox(ue_flows, eq_flows)
      end

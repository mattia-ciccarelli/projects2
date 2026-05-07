
# inflow condition fluidic pinball
# function inflow(coor::Vector{Float64},idof::Int64)

#   if idof==2
#     valdbc=0.0
#   else
#     valdbc=Uin
#   end

#   return valdbc
# end

# inflow condition flow around cylinder
function inflow(coor::Vector{Float64},idof::Int64)

  # H=0.41
  # H=0.4
  # H=4.1
  if idof==2
    valdbc=0.0
  else
    valdbc=Uin
    # valdbc=4*Uin/H^2*coor[2]*(H-coor[2])
  end

return valdbc
end



if !("FPBH" in keys(Pkg.installed()))
	Pkg.clone("https://github.com/aritrasep/FPBH.jl")
	Pkg.build("FPBH")
end

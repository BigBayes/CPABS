macro ifndef(sym, val)
    if !isdefined(sym)
        return esc(:($sym = $val))
    end
end

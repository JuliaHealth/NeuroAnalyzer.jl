export na_test_plugin

"""
    na_test_plugin()

NeuroAnalyzer test plugin.

# Arguments

Nothing

# Returns

Nothing
"""
function na_test_plugin()::Nothing
    @info "This is a test plugin for NeuroAnalyzer"
    @info "Running na_info()"
    na_info()
    @info "Test completed"

    return nothing

end

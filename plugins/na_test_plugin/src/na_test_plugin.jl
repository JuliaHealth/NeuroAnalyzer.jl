export na_test_plugin

"""
    na_test_plugin()

NeuroAnalyzer test plugin.

# Arguments

Nothing

# Returns

- `Nothing`
"""
function na_test_plugin()::Nothing

    _info("This is a test plugin for NeuroAnalyzer")
    _info("Running na_info()")

    na_info()

    _info("Test completed")

    return nothing

end

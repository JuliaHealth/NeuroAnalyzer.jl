export test_plugin

"""
    test_plugin()

NeuroAnalyzer test plugin.

# Arguments

Nothing

# Returns

- `Nothing`
"""
function test_plugin()::Nothing
    _info("This is a test plugin for NeuroAnalyzer")
    _info("Running info()")

    info()

    _info("Test completed")

    return nothing
end

export na_test_plugin

"""
    na_test_plugin()

NeuroAnalyzer test plugin.
"""
function na_test_plugin()
    @info "This is a test plugin for NeuroAnalyzer"
    @info "Running na_info()"
    na_info()
    @info "Test completed"
end
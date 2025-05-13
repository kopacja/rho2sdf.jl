using Dates

# Structure to hold logger configuration
struct LoggerConfig
    log_file::String
    log_level::Symbol
    include_timestamp::Bool
    write_to_console::Bool
end

# ANSI color codes for different message types
const COLORS = Dict(
    :error => "\e[1;31m",    # Red
    :warning => "\e[1;33m",  # Yellow
    :info => "\e[1;34m",     # Blue
    :success => "\e[1;32m"   # Green
)

# Create a new log file or append to existing one
function initialize_logger(filename::String="computation_log.txt",
                         log_level::Symbol=:info,
                         include_timestamp::Bool=true,
                         write_to_console::Bool=true)
    
    # Create the logs directory if it doesn't exist
    log_dir = "logs"
    !isdir(log_dir) && mkdir(log_dir)
    
    # Add timestamp to filename to make it unique
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    full_filename = joinpath(log_dir, "$(timestamp)_$(filename)")
    
    # Create the log file with a header
    open(full_filename, "a") do file
        println(file, "=== Computation Log Started at $(timestamp) ===")
        println(file, "Log Level: $(log_level)")
        println(file, "================================================")
    end
    
    return LoggerConfig(full_filename, log_level, include_timestamp, write_to_console)
end

# Main logging function
function log_message(config::LoggerConfig, level::Symbol, message::String)
    # Define log level hierarchy
    level_hierarchy = Dict(
        :error => 1,
        :warning => 2,
        :info => 3,
        :success => 4
    )
    
    # Check if we should log this message based on log level
    if level_hierarchy[level] > level_hierarchy[config.log_level]
        return
    end
    
    # Create timestamp if needed
    timestamp = config.include_timestamp ? "[$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))] " : ""
    
    # Format the message
    formatted_msg = "$(timestamp)[$(uppercase(String(level)))] $(message)"
    
    # Write to log file
    open(config.log_file, "a") do file
        println(file, formatted_msg)
    end
    
    # Write to console if enabled
    if config.write_to_console
        color = get(COLORS, level, "")
        println(color, formatted_msg, "\e[0m")
    end
end

# Convenience functions for different log levels
function log_error(config::LoggerConfig, message::String)
    log_message(config, :error, message)
end

function log_warning(config::LoggerConfig, message::String)
    log_message(config, :warning, message)
end

function log_info(config::LoggerConfig, message::String)
    log_message(config, :info, message)
end

function log_success(config::LoggerConfig, message::String)
    log_message(config, :success, message)
end

# Function to log computation progress
function log_computation_progress(config::LoggerConfig, 
                                element_number::Int, 
                                total_elements::Int, 
                                current_volume::Float64)
    progress = round(element_number / total_elements * 100, digits=2)
    log_info(config, "Progress: $(progress)% - Element: $(element_number)/$(total_elements) - Current Volume: $(current_volume)")
end


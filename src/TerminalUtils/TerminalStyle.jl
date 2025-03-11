# Info messages - blue color for general information
function print_info(msg::String)
    # Count leading newlines
    leading_newlines = length(match(r"^\n*", msg).match)
    # Remove leading newlines and get clean message
    cleaned_msg = lstrip(msg, '\n')
    
    # Print the leading newlines first if any exist
    if leading_newlines > 0
        print("\n" ^ leading_newlines)
    end
    
    # Print the formatted message
    println("\e[1;34m[INFO] $(cleaned_msg)\e[0m")
end

# Error messages - red color for critical issues
function print_error(msg::String)
    leading_newlines = length(match(r"^\n*", msg).match)
    cleaned_msg = lstrip(msg, '\n')
    
    if leading_newlines > 0
        print("\n" ^ leading_newlines)
    end
    println("\e[1;31m[ERROR] $(cleaned_msg)\e[0m")
end

# Warning messages - orange/yellow color for potential issues
function print_warning(msg::String)
    leading_newlines = length(match(r"^\n*", msg).match)
    cleaned_msg = lstrip(msg, '\n')
    
    if leading_newlines > 0
        print("\n" ^ leading_newlines)
    end
    println("\e[1;33m[WARNING] $(cleaned_msg)\e[0m")
end

# Success messages - green color for completion and positive results
function print_success(msg::String)
    leading_newlines = length(match(r"^\n*", msg).match)
    cleaned_msg = lstrip(msg, '\n')
    
    if leading_newlines > 0
        print("\n" ^ leading_newlines)
    end
    println("\e[1;32m[SUCCESS] $(cleaned_msg)\e[0m")
end

# Highlighting important messages
function print_data(msg::String)
    leading_newlines = length(match(r"^\n*", msg).match)
    cleaned_msg = lstrip(msg, '\n')
    
    if leading_newlines > 0
        print("\n" ^ leading_newlines)
    end
    println("\e[33m $(cleaned_msg)\e[0m")
end

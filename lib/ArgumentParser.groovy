/* Check arguments of a Nextflow function
 *
 * Nextflow script does not support the Groovy idiom:
 *
 *     def function(Map args[:], arg1, arg2, ...)
 * 
 * to support unordered kwargs. The methods here are designed
 * to reduce boileplate while allowing Nextflow script to implement
 *
 *     def function(Map args[:])
 *
 * with required and default values. This is similar to some Python
 * libraries' (notably matplotlib) extensive use of things like:
 *
 *     def function(*args, **kwargs)
 *
 * to implement generic APIs. Why do we want to do all this? Because
 * we want to write library code with a clean set of required parameters
 * but also extensible with non-required parameters with default values.
 * This allows us to later add parameters without breaking existing code,
 * and is very common practice elsewhere.
 */

import java.util.Set

class ArgumentParser {
    Set args
    Map kwargs
    String name

    /* Parse arguments, raising an error on unknown keys */
    public Map parse_args(LinkedHashMap given_args) {
        Set opt_keys = kwargs.keySet()
        Set given_keys = given_args.keySet()
        check_required(given_keys)
        check_unknown(given_keys, opt_keys)
        return kwargs + given_args
    }
    
    /* Parse arguments, without raising an error for extra keys */
    public Map parse_known_args(LinkedHashMap given_args) {
        Set opt_keys = kwargs.keySet()
        Set given_keys = given_args.keySet()
        check_required(given_keys)
        return kwargs + given_args
    }
    
    private void check_required(Set given) {
        Set missing_keys = args - given
        if (!missing_keys.isEmpty()) {
            throw new Exception("Missing arguments for function ${name}: ${missing_keys}")
        }
    }
    
    private void check_unknown(Set given, Set kwargs_keys) {
        Set extra_keys = given - (args + kwargs_keys)
        if (!extra_keys.isEmpty()) {
            throw new Exception("Unknown arguments provided to function ${name}: ${extra_keys}.")
        }
    }
}

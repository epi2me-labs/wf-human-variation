/* Miscellaneous utilities for workflows from the ONT Customer Workflows Group.
 */
class CWUtil {

    /* Mutate the global Nextflow params map
    *
    * Occasionally, we may wish to mutate the value of a parameter provided
    * by the user. Typically, this leads to workflows with `params.my_param`
    * and `params._my_param` which is ripe for confusion. Instead, we can
    * mutate the parameter value in the Nextflow params ScriptMap itself
    * with the following call:
    *
    *     CWUtil.mutateParam(params, k, v)
    *
    * This is possible as Groovy actually has a surprisingly loose
    * definition of "private", and allows us to call the private `allowNames`
    * method on the ScriptMap which removes the read-only status for a key set.
    * We can follow this up with a call to the private `put0` to reinsert
    * the key and mark it as read-only again.
    */
    public static void mutateParam(nf_params, key, value) {
        Set s = [key] // must be a set to allow call to allowNames
        nf_params.allowNames(s)
        nf_params.put0(key, value)
    }
}

import random
import six

# to calculate the standard Robinson-Foulds (RF) distance
def robinson_foulds(t1, t2, attr_t1="name", attr_t2="name"):
    '''compute Robinson Foulds Distance,
       supporting non-binary tree

       Parameters
       ----------
       t1, t2: ete3 Tree like objects

       Returns
       -------
       results: dic objects
       results['rf'] Robinson Foulds Distance
       results['norm_rf] Normalized Robinson Foulds Distance
    '''
    ref_t = t1
    target_t = t2



    attrs_t1 = set([getattr(n, attr_t1) for n in ref_t.iter_leaves() if hasattr(n, attr_t1)])
    attrs_t2 = set([getattr(n, attr_t2) for n in target_t.iter_leaves() if hasattr(n, attr_t2)])
    common_attrs = attrs_t1 & attrs_t2
    # release mem
    attrs_t1, attrs_t2 = None, None

    # Check for duplicated items (is it necessary? can we optimize? what's the impact in performance?')
    size1 = len([True for n in ref_t.iter_leaves() if getattr(n, attr_t1, None) in common_attrs])
    size2 = len([True for n in target_t.iter_leaves() if getattr(n, attr_t2, None) in common_attrs])
    if size1 > len(common_attrs):
        raise TreeError('Duplicated items found in source tree')
    if size2 > len(common_attrs):
        raise TreeError('Duplicated items found in reference tree')

    t1_content = t1.get_cached_content()
      
    edges1 = set([
            tuple(sorted([getattr(n, attr_t1) for n in content if hasattr(n, attr_t1) and getattr(n, attr_t1) in common_attrs]))
            for content in six.itervalues(t1_content)])
    edges1.discard(())
    
    t2_content = t2.get_cached_content()

    edges2 = set([
            tuple(sorted([getattr(n, attr_t2) for n in content if hasattr(n, attr_t2) and getattr(n, attr_t2) in common_attrs]))
            for content in six.itervalues(t2_content)])
    edges2.discard(())

    # the two root edges are never counted here, as they are always
    # present in both trees because of the common attr filters
    rf = len((edges1 ^ edges2))

    max_parts = len([p for p in edges1 if len(p)>1]) 
    max_parts +=  len([p for p in edges2 if len(p)>1]) - 2

    result = {}
    result['rf'] = rf
    result["rf_norm"] = rf / max_parts

    return result 
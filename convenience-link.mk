# Makes links to the hidden libtool-built libraries in the build directory as
# normal non-libtool builds usually make.
#
# Derived from a StackOverflow comment: http://stackoverflow.com/a/14219698
#

.PHONY: convenience-link clean-convenience-link

convenience-link: $(lib_LTLIBRARIES) $(noinst_LTLIBRARIES)
	$(AM_V_at)for soname in `echo | $(EGREP) "^(dlname|old_library)=" $^ | $(SED) -e "s#^\(dlname\|old_library\)='\(.*\)'#\2#"`; do  \
		rm -f $(builddir)/$$soname; \
		test -e $(builddir)/.libs/$$soname && \
		cd $(builddir) && \
		$(LN_S) $(builddir)/.libs/$$soname $$soname || true;\
	done

clean-convenience-link:
	$(AM_V_at)for linkname in `ls -1 $(builddir) | $(EGREP) '^$(lib_LTLIBRARIES:.la=.)'`; do \
		echo "test -L $(builddir)/$$linkname && rm -f $(builddir)/$$linkname"; \
		test -L $(builddir)/$$linkname && rm -f $(builddir)/$$linkname || true; \
	done

all-local:: convenience-link

clean-local:: clean-convenience-link


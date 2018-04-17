package org.jcvi.vigor.service;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.component.MaturePeptideMatch;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.VigorUtils;
import org.junit.ClassRule;
import org.junit.Rule;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.rules.SpringClassRule;
import org.springframework.test.context.junit4.rules.SpringMethodRule;

import java.io.File;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.CoreMatchers.hasItem;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;

@RunWith(Parameterized.class)
@ContextConfiguration(classes = Application.class)
public class PeptideServiceTest {

    private final static Logger LOGGER = LogManager.getLogger(PeptideServiceTest.class);

    @ClassRule
    public static final SpringClassRule springClassRule = new SpringClassRule();

    @Rule
    public final SpringMethodRule springMethodRule = new SpringMethodRule();


    @Autowired
    private PeptideService peptideService;

    private Collection<String> expected;

    private String id;
    private String sequence;
    private String mp_ref_db;

    public PeptideServiceTest(String id, String sequence, String mp_ref_db, Collection<String> expected) {
        this.id = id;
        this.sequence = sequence;
        this.mp_ref_db = mp_ref_db;
        this.expected = expected;
    }

    @Test
    public void testPeptides() throws ServiceException {

        // must have atleast id, polyprotein, peptide db
        ViralProtein protein = new ViralProtein();
        protein.setDefline(id);
        protein.setProteinID("");
        protein.setSequence(ProteinSequence.of(sequence));
        File peptideDB = new File(mp_ref_db);

        List<MaturePeptideMatch> matches = peptideService.findPeptides(protein, peptideDB);
        LOGGER.debug(() -> String.format("peptides:%s", matches.stream().map(String::valueOf).collect(Collectors.joining("\n> ", "\n> ", ""))));

        assertThat(String.format("peptide mismatches for %s", protein.getDefline()), matches.size(), equalTo(expected.size()));
        List<ProteinSequence> subjectMatches = matches.stream()
                                                      .map(m -> m.getProtein().getSequence().toBuilder().trim(m.getProteinRange()).build())
                                                      .collect(Collectors.toList());

        for (String expectedPeptide : expected) {
            assertThat(subjectMatches, hasItem(ProteinSequence.of(expectedPeptide)));
        }
        Range prev = null;
        Range current;
        MaturePeptideMatch prevMatch = null;
        for (MaturePeptideMatch match: matches) {
            current = match.getProteinRange();
            if (prev == null) {
                assertTrue(current.getBegin() == 0 || match.isFuzzyBegin());
            } else {
                assertTrue (String.format("ranges [%s:%s] fuzzyEnd=%s, fuzzyBegin=%s [%s:%s]  are not contiguous and not fuzzy",
                        prev.getBegin(), prev.getEnd(), prevMatch.isFuzzyEnd(),
                        match.isFuzzyBegin(), current.getBegin(), current.getEnd()), current.getBegin() == prev.getEnd() + 1 ||
                        ( match.isFuzzyBegin() && prevMatch.isFuzzyEnd() ));
            }
            prev = current;
            prevMatch = match;
        }
        // TODO test produced names and products

    }


    @Parameterized.Parameters(name="PeptideServiceTest[#{index} {0}]")
    public static Collection<Object[]> getTestData() {

        List<Object[]> testData = new ArrayList<>();
        testData.add(
                new Object[]{
                        "HA_1_615",
                        String.join("",
                                "MKAILVVLLYTFATANADTLCIGYHANNST",
                                "DTVDTVLEKNVTVTHSVNLLEDKHNGKLCK",
                                "LRGVAPLHLGKCNIAGWILGNPECESLSTA",
                                "SSWSYIVETSSSDNGTCYPGDFINYEELRE",
                                "QLSSVSSFERFEIFPKTSSWPNHDSNKGVT",
                                "AACPQAGAKSFYKNLIWLVKKGNSYPKLSK",
                                "SYINDKGKEVLVLWGIHH"),
                        "flua_ha_mp",
                        Arrays.asList(
                        "MKAILVVLLYTFATANA",
                        String.join("",
                                "DTLCIGYHANNSTDTVDTVLEKNVTVTHSV",
                                "NLLEDKHNGKLCKLRGVAPLHLGKCNIAGW",
                                "ILGNPECESLSTASSWSYIVETSSSDNGTC",
                                "YPGDFINYEELREQLSSVSSFERFEIFPKT",
                                "SSWPNHDSNKGVTAACPQAGAKSFYKNLIW",
                                "LVKKGNSYPKLSKSYINDKGKEVLVLWGIHH"))
                });
        testData.add(
                new Object[]{
                        "HA_1129_1753.1",
                        String.join("",
                                "QNAQGSGYAADLKSTQNAIDKITNKVNSVI",
                                "EKMNTQFTAVGKEFNHLEKRIENLNKKVDD",
                                "GFLDIWTYNAELLVLLENERTLDYHDSNVK",
                                "NLYEKVRSQLKNNAKEIGNGCFEFYHKCDN",
                                "TCMESVKNGTYDYPRYSEEAKLNREEIDGV",
                                "KLESTRIYQILAIYSTVASSLVLIVSLGAI",
                                "SFWMCSNGSLQCRICI*"),
                        "flua_ha_mp",
                        Arrays.asList(
                        String.join("",
                                "QNAQGSGYAADLKSTQNAIDKITNKVNSVI",
                                "EKMNTQFTAVGKEFNHLEKRIENLNKKVDD",
                                "GFLDIWTYNAELLVLLENERTLDYHDSNVK",
                                "NLYEKVRSQLKNNAKEIGNGCFEFYHKCDN",
                                "TCMESVKNGTYDYPRYSEEAKLNREEIDGV",
                                "KLESTRIYQILAIYSTVASSLVLIVSLGAI",
                                "SFWMCSNGSLQCRICI"))
                });
        testData.add(
                new Object[]{
                        "gi_260907762.1",
                        String.join("",
                                "MNIQILVFALVAIIPTNADKICLGHHAVSN",
                                "GTKVNTLTERGVEVVNATETVERTNVPRIC",
                                "SKGKRTVDLGQCGLLGTITGPPQCDQFLEF",
                                "SADLIIERRGGSDVCYPGKFVNEEALRQIL",
                                "RESGGIDKETMGFTYSGIRTNGATSACRRS",
                                "GSSFYAEMKWLLSNTDNAAFPQMTKSYKNT",
                                "RKDPALIIWGIHHSGSTTEQTKLYGSGSKL",
                                "ITVGSSNYQQSFVPSPGARPQVNGQSGRID",
                                "FHWLMLNPNDTVTFSFNGAFIAPDRASFLK",
                                "GKSMGIQSGVQVDANCEGDCYHSGGTIISN",
                                "LPFQNINSRAVGKCPRYVKQESLMLATGMK",
                                "NVPELPKGRGLFGAIAGFIENGWEGLIDGW",
                                "YGFRHQNAQGEGTAADYKSTQSAIDQITGK",
                                "LNRLIEKTNQQFELIDNEFTEVEKQIGNVI",
                                "NWTRDSMTEVWSYNAELLVAMENQHTIDLA",
                                "DSEMNKLYERVKRQLRENAEEDGTGCFEIF",
                                "HKCDDDCMASIRNNTYDHSKYREEAMQNRI",
                                "QINPVKLSSGYKDVILWFSFGASCFILLAI",
                                "AMGLVFICVKNGNMRCTICI*"),
                        "flua_ha_mp",
                        Arrays.asList(
                        "MNIQILVFALVAIIPTNA",
                        String.join("",
                                "DKICLGHHAVSNGTKVNTLTERGVEVVNAT",
                                "ETVERTNVPRICSKGKRTVDLGQCGLLGTI",
                                "TGPPQCDQFLEFSADLIIERRGGSDVCYPG",
                                "KFVNEEALRQILRESGGIDKETMGFTYSGI",
                                "RTNGATSACRRSGSSFYAEMKWLLSNTDNA",
                                "AFPQMTKSYKNTRKDPALIIWGIHHSGSTT",
                                "EQTKLYGSGSKLITVGSSNYQQSFVPSPGA",
                                "RPQVNGQSGRIDFHWLMLNPNDTVTFSFNG",
                                "AFIAPDRASFLKGKSMGIQSGVQVDANCEG",
                                "DCYHSGGTIISNLPFQNINSRAVGKCPRYV",
                                "KQESLMLATGMKNVPELPKGR"),
                        String.join("",
                                "GLFGAIAGFIENGWEGLIDGWYGFRHQNAQ",
                                "GEGTAADYKSTQSAIDQITGKLNRLIEKTN",
                                "QQFELIDNEFTEVEKQIGNVINWTRDSMTE",
                                "VWSYNAELLVAMENQHTIDLADSEMNKLYE",
                                "RVKRQLRENAEEDGTGCFEIFHKCDDDCMA",
                                "SIRNNTYDHSKYREEAMQNRIQINPVKLSS",
                                "GYKDVILWFSFGASCFILLAIAMGLVFICV",
                                "KNGNMRCTICI"))
                });
        testData.add(
                new Object[]{
                        "gi|260907760|gb|GU060481.1|",
                        String.join("",
                                "MNTQILALIACMLIGAKGDKICLGHHAVAN",
                                "GTKVNTLTERGIEVVNATETVETANIKKIC",
                                "TQGKRPTDLGQCGLLGTLIGPPQCDQFLEF",
                                "DTDLIIERREGTDVCYPGKFTNEESLRQIL",
                                "RGSGGIDKESMGFTYSGIRTNGATSACRRS",
                                "GSSFYAEMKWLLSNSDNAAFPQMTKSYRNP",
                                "RNKPALIIWGVHHSGSATEQTKLYGSGNKL",
                                "ITVGSSKYQQSFTPSPGARPQVNGQSGRID",
                                "FHWLLLDPNDTVTFTFNGAFIAPDRASFFR",
                                "GESLGVQSDVPLDSGCEGDCFHSGGTIVSS",
                                "LPFQNINPRTVGKCPRYVKQTSLLLATGMR",
                                "NVPENPKTRGLFGAIAGFIENGWEGLIDGW",
                                "YGFRHQNAQGEGTAADYKSTQSAIDQITGK",
                                "LNRLIDKTNQQFELIDNEFSEIEQQIGNVI",
                                "NWTRDSMTEVWSYNAELLVAMENQHTIDLA",
                                "DSEMNKLYERVRKQLRENAEEDGTGCFEIF",
                                "HKCDDQCMESIRNNTYDHTQYRTESLQNRI",
                                "QIDPVRLSSGYKDIILWFSFGASCFLLLAI",
                                "AMGLVFICIKNGNMRCTICI*"),
                        "flua_ha_mp",
                        Arrays.asList(
                        "MNTQILALIACMLIGAKG",
                        String.join("",
                                "DKICLGHHAVANGTKVNTLTERGIEVVNAT",
                                "ETVETANIKKICTQGKRPTDLGQCGLLGTL",
                                "IGPPQCDQFLEFDTDLIIERREGTDVCYPG",
                                "KFTNEESLRQILRGSGGIDKESMGFTYSGI",
                                "RTNGATSACRRSGSSFYAEMKWLLSNSDNA",
                                "AFPQMTKSYRNPRNKPALIIWGVHHSGSAT",
                                "EQTKLYGSGNKLITVGSSKYQQSFTPSPGA",
                                "RPQVNGQSGRIDFHWLLLDPNDTVTFTFNG",
                                "AFIAPDRASFFRGESLGVQSDVPLDSGCEG",
                                "DCFHSGGTIVSSLPFQNINPRTVGKCPRYV",
                                "KQTSLLLATGMRNVPENPKT"),
                        String.join("",
                                "GLFGAIAGFIENGWEGLIDGWYGFRHQNAQ",
                                "GEGTAADYKSTQSAIDQITGKLNRLIDKTN",
                                "QQFELIDNEFSEIEQQIGNVINWTRDSMTE",
                                "VWSYNAELLVAMENQHTIDLADSEMNKLYE",
                                "RVRKQLRENAEEDGTGCFEIFHKCDDQCME",
                                "SIRNNTYDHTQYRTESLQNRIQIDPVRLSS",
                                "GYKDIILWFSFGASCFLLLAIAMGLVFICI",
                                "KNGNMRCTICI"
                        ))
                });
        return testData;
    }


}

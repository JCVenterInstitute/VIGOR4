package org.jcvi.vigor.service;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.component.MaturePeptideMatch;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.VigorConfiguration;
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
import java.net.URL;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.CoreMatchers.hasItem;
import static org.hamcrest.CoreMatchers.notNullValue;
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
    @Autowired
    private VigorInitializationService initializationService;

    private List<String[]> expected;

    private String id;
    private String sequence;
    private String mp_ref_db;

    public PeptideServiceTest(String id, String sequence, String mp_ref_db, List<String[]> expected) {
        this.id = id;
        this.sequence = sequence;
        this.mp_ref_db = mp_ref_db;
        this.expected = expected;
    }

    @Test
    public void testPeptides() throws VigorException {

        // must have atleast id, polyprotein, peptide db
        ProteinSequence protein = ProteinSequence.of(sequence);
        VigorConfiguration config = initializationService.mergeConfigurations(initializationService.getDefaultConfigurations());
        String refDBPath = config.get(ConfigurationParameters.ReferenceDatabasePath);
        assertThat("Reference database path must be set", refDBPath, notNullValue());
        String refDB = Paths.get(refDBPath, mp_ref_db).toString();
        File peptideDB = getPeptideDB(refDBPath, mp_ref_db);

        List<MaturePeptideMatch> matches = peptideService.findPeptides(protein, peptideDB);
        LOGGER.debug(() -> String.format("peptides:%s", matches.stream().map(String::valueOf).collect(Collectors.joining("\n> ", "\n> ", ""))));

        assertThat(String.format("peptide mismatches for %s", id), matches.size(), equalTo(expected.size()));

        assertThat("matches count should be that same as expected count",matches.size(), equalTo(expected.size()));

        String expectedProduct;
        String expectedID;
        ProteinSequence expectedSequence;
        MaturePeptideMatch match;
        ProteinSequence matchedSequence;
        Range prev = null;
        Range current;
        MaturePeptideMatch prevMatch = null;

        for (int i = 0; i < expected.size(); i++) {
            match = matches.get(i);
            expectedID = expected.get(i)[0];
            expectedProduct = expected.get(i)[1];
            expectedSequence = ProteinSequence.of(expected.get(i)[2]);
            matchedSequence =  match.getProtein().toBuilder().trim(match.getProteinRange()).build();
            assertThat(matchedSequence, equalTo(expectedSequence));
            assertThat(match.getReference().getProduct(), equalTo(expectedProduct));

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

    }

    private File getPeptideDB(String referenceDBPath, String mp_ref_db) throws VigorException {
        File peptideDB = Paths.get(referenceDBPath,mp_ref_db).toFile();
        if (! peptideDB.exists()) {
            throw new VigorException(String.format("unable to find peptide DB for %s", mp_ref_db));
        }
        return peptideDB;
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
                                new String[] {"seg4matureA1","signal peptide", "MKAILVVLLYTFATANA" },
                                new String[] { "seg4matureA2","HA1",
                        String.join("",
                                "DTLCIGYHANNSTDTVDTVLEKNVTVTHSV",
                                "NLLEDKHNGKLCKLRGVAPLHLGKCNIAGW",
                                "ILGNPECESLSTASSWSYIVETSSSDNGTC",
                                "YPGDFINYEELREQLSSVSSFERFEIFPKT",
                                "SSWPNHDSNKGVTAACPQAGAKSFYKNLIW",
                                "LVKKGNSYPKLSKSYINDKGKEVLVLWGIHH") })
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
                        Arrays.<String[]>asList( new String[] {"seg4matureA3","HA2",
                        String.join("",
                                "QNAQGSGYAADLKSTQNAIDKITNKVNSVI",
                                "EKMNTQFTAVGKEFNHLEKRIENLNKKVDD",
                                "GFLDIWTYNAELLVLLENERTLDYHDSNVK",
                                "NLYEKVRSQLKNNAKEIGNGCFEFYHKCDN",
                                "TCMESVKNGTYDYPRYSEEAKLNREEIDGV",
                                "KLESTRIYQILAIYSTVASSLVLIVSLGAI",
                                "SFWMCSNGSLQCRICI") })
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
                        Arrays.asList( new String[] { "seg4matureO1","signal peptide",
                        "MNIQILVFALVAIIPTNA" },
                                new String[] { "seg4mature7F2", "HA1",
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
                                "KQESLMLATGMKNVPELPKGR") },
                                new String[] { "seg4matureO3","HA2",
                        String.join("",
                                "GLFGAIAGFIENGWEGLIDGWYGFRHQNAQ",
                                "GEGTAADYKSTQSAIDQITGKLNRLIEKTN",
                                "QQFELIDNEFTEVEKQIGNVINWTRDSMTE",
                                "VWSYNAELLVAMENQHTIDLADSEMNKLYE",
                                "RVKRQLRENAEEDGTGCFEIFHKCDDDCMA",
                                "SIRNNTYDHSKYREEAMQNRIQINPVKLSS",
                                "GYKDVILWFSFGASCFILLAIAMGLVFICV",
                                "KNGNMRCTICI") })
                });
        testData.add(
                new Object[]{
                        "gi_155016323.1",
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
                        Arrays.asList( new String[] { "seg4mature7A1","signal peptide",
                        "MNTQILALIACMLIGAKG" },
                                new String[] { "seg4mature7C2","HA1",
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
                                "KQTSLLLATGMRNVPENPKTR") },
                                new String[] { "seg4mature7A3", "HA2",
                        String.join("",
                                "GLFGAIAGFIENGWEGLIDGWYGFRHQNAQ",
                                "GEGTAADYKSTQSAIDQITGKLNRLIDKTN",
                                "QQFELIDNEFSEIEQQIGNVINWTRDSMTE",
                                "VWSYNAELLVAMENQHTIDLADSEMNKLYE",
                                "RVRKQLRENAEEDGTGCFEIFHKCDDQCME",
                                "SIRNNTYDHTQYRTESLQNRIQIDPVRLSS",
                                "GYKDIILWFSFGASCFLLLAIAMGLVFICI",
                                "KNGNMRCTICI") }
                        )
                });
        return testData;
    }


}

from sqlalchemy import (
    Boolean,
    CheckConstraint,
    Column,
    Date,
    DateTime,
    Float,
    ForeignKey,
    Integer,
    String,
    Table,
    event,
    text,
    true,
    orm,
    select,
    update,
)
import sqlalchemy.sql as sql
from sqlalchemy.sql.sqltypes import NullType
from sqlalchemy.orm import Session, relationship, backref
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm.collections import column_mapped_collection
from itertools import chain
import datetime

NEXT_ID_SUBQUERY = "SELECT coalesce(max(id)+1, 0) FROM"

Base = declarative_base()
metadata = Base.metadata

# TODO: event listeners get duplicated when this file is autoreloaded, this mangles queries

# FROM: https://docs.sqlalchemy.org/en/20/_modules/examples/extending_query/filter_public.html
@event.listens_for(Session, "do_orm_execute")
def only_show_visible(execute_state):
    """Intercept all ORM queries.   Add a with_loader_criteria option to all of
    them.

    This option applies to SELECT queries and adds a global WHERE
    criteria (or as appropriate ON CLAUSE criteria for join targets) to
    all objects of a certain class or superclass.
    """

    # the with_loader_criteria automatically applies itself to
    # relationship loads as well including lazy loads.   So if this is
    # a relationship load, assume the option was set up from the top level
    # query.

    if (
        not execute_state.is_column_load
        and not execute_state.is_relationship_load
        and not execute_state.execution_options.get("include_invisible", False)
    ):
        execute_state.statement = execute_state.statement.options(
            orm.with_loader_criteria(
                HasVisible, lambda cls: cls.visible == true(), include_aliases=True
            )
        )


class HasVisible:
    visible = Column(Boolean, nullable=False, default=True)


@event.listens_for(Session, "after_flush")
def update_next_table_id(session, flush_context):
    klasses = set(
        obj.__class__
        for obj in chain(session.new, session.dirty)
        if hasattr(obj.__class__, "id")
    )
    # for table in tables:
    for klass in klasses:
        next_id_subquery = select(sql.functions.max(klass.id) + 1).scalar_subquery()
        session.execute(
            update(NextTableId)
            .where(NextTableId.table_name == klass.__tablename__)
            .values(next_id=next_id_subquery)
        )
    #     session.execute(sql.functions.max(g.Table.id))
    # print(tables, dir(session), "||", dir(flush_context))
    # session.execute("UPDATE ")


class DocumentFileDatum(Base):
    __tablename__ = "document_file_data"

    id = Column(
        Integer, primary_key=True, default=text(f"({NEXT_ID_SUBQUERY} {__tablename__})")
    )
    data = Column(NullType)
    local_file_path = Column(String)
    local_file_size = Column(Integer)
    last_needed = Column(DateTime)


class GGroup(Base):
    __tablename__ = "g_group"

    id = Column(
        Integer, primary_key=True, default=text(f"({NEXT_ID_SUBQUERY} {__tablename__})")
    )
    name = Column(String(255), nullable=False)

    def __repr__(self):
        return f"<GGroup '{self.name}' id={self.id}>"


class GRole(Base):
    __tablename__ = "g_role"

    id = Column(
        Integer, primary_key=True, default=text(f"({NEXT_ID_SUBQUERY} {__tablename__})")
    )
    name = Column(String(255), nullable=False)

    def __repr__(self):
        return f"<GRole '{self.name}' id={self.id}>"


class Metadata(Base):
    __tablename__ = "metadata"

    identifier = Column(String(80), primary_key=True)
    value = Column(String(255))

    def __repr__(self):
        return f"<Metadata '{self.identifier}'='{self.value}'>"


class NextTableId(Base):
    __tablename__ = "next_table_id"

    table_name = Column(String(50), primary_key=True)
    next_id = Column(Integer)

    def __repr__(self):
        return f"<NextTableId {self.table_name}: {self.next_id}>"


class SearchField(Base):
    __tablename__ = "search_field"

    code = Column(String(255), primary_key=True)
    field_xml = Column(String, nullable=False)


class Folder(HasVisible, Base):
    __tablename__ = "folder"
    __table_args__ = (CheckConstraint("id != parent_folder_id"),)

    id = Column(
        Integer, primary_key=True, default=text(f"({NEXT_ID_SUBQUERY} {__tablename__})")
    )
    g_group_id = Column(ForeignKey("g_group.id"), nullable=False)
    parent_folder_id = Column(ForeignKey("folder.id"), index=True)
    # modified = Column(DateTime, nullable=False, default=datetime.datetime.utcnow)
    # modified = Column(DateTime, nullable=False, default=datetime.datetime.now)
    modified = Column(DateTime, nullable=False, default=text("current_timestamp"))
    name = Column(String(255), index=True)

    g_group = relationship("GGroup")
    children = relationship(
        "Folder",
        cascade="all, delete-orphan",
        backref=backref("parent_folder", remote_side=[id]),
        collection_class=column_mapped_collection(name),
    )
    annotated_documents = relationship(
        "AnnotatedDocument",
        cascade="all, delete-orphan",
        backref=backref("parent_folder"),
    )
    users = relationship("GUser", secondary="hidden_folder_to_user")

    def __contains__(self, key):
        return key in self.children

    def __getitem__(self, key):
        return self.children[key]

    def __repr__(self):
        return f"<Folder '{self.name}' id={self.id} parent_folder_id={self.parent_folder_id}>"


class GUser(Base):
    __tablename__ = "g_user"

    id = Column(
        Integer, primary_key=True, default=text(f"({NEXT_ID_SUBQUERY} {__tablename__})")
    )
    primary_group_id = Column(ForeignKey("g_group.id"), nullable=False, index=True)
    username = Column(String(255), nullable=False)

    primary_group = relationship("GGroup")

    def __repr__(self):
        return f"<GUser '{self.username}' id={self.id} primary_group_id={self.primary_group_id}>"


class AnnotatedDocument(Base):
    __tablename__ = "annotated_document"

    id = Column(
        Integer, primary_key=True, default=text(f"({NEXT_ID_SUBQUERY} {__tablename__})")
    )
    folder_id = Column(
        ForeignKey("folder.id", ondelete="CASCADE"), nullable=False, index=True
    )
    modified = Column(DateTime, nullable=False)
    urn = Column(String(255), nullable=False, unique=True)
    document_xml = Column(String, nullable=False)
    plugin_document_xml = Column(String, nullable=False)
    reference_count = Column(Integer, nullable=False)

    # folder = relationship("Folder")
    g_users = relationship("GUser", secondary="document_read")
    file_datas = relationship("DocumentFileDatum", secondary="document_to_file_data")
    additional_document_xml = relationship(
        "AdditionalDocumentXml",
        cascade="all, delete-orphan",
        backref=backref("annotated_document"),
    )
    additional_xml_timestamp = relationship(
        "AdditionalXmlTimestamp",
        cascade="all, delete-orphan",
        backref=backref("annotated_document"),
        uselist=False,
    )

    def __repr__(self):
        return f"<AnnotatedDocument id={self.id} {self.urn} modified={self.modified}>"


class IndexingQueue(AnnotatedDocument):
    __tablename__ = "indexing_queue"

    document_id = Column(
        ForeignKey("annotated_document.id", ondelete="CASCADE"), primary_key=True
    )
    g_user_id = Column(ForeignKey("g_user.id", ondelete="SET NULL"), index=True)
    reserved = Column(DateTime)

    g_user = relationship("GUser")


class FolderView(Base):
    __tablename__ = "folder_view"

    folder_id = Column(
        ForeignKey("folder.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    )
    document_urn = Column(String(255), primary_key=True, nullable=False)
    modified = Column(DateTime, nullable=False)

    folder = relationship("Folder")


class GUserGroupRole(Base):
    __tablename__ = "g_user_group_role"

    g_user_id = Column(
        ForeignKey("g_user.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    )
    g_group_id = Column(
        ForeignKey("g_group.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    )
    g_role_id = Column(
        ForeignKey("g_role.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    )

    g_group = relationship("GGroup")
    g_role = relationship("GRole")
    g_user = relationship("GUser")


t_hidden_folder_to_user = Table(
    "hidden_folder_to_user",
    metadata,
    Column(
        "hidden_folder_id",
        ForeignKey("folder.id", ondelete="CASCADE"),
        primary_key=True,
    ),
    Column("user_id", ForeignKey("g_user.id", ondelete="CASCADE"), index=True),
)


class SpecialElement(Base):
    __tablename__ = "special_element"

    folder_id = Column(
        ForeignKey("folder.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    )
    modified = Column(DateTime, nullable=False)
    xml = Column(String, nullable=False)
    name = Column(String(255), primary_key=True, nullable=False)

    folder = relationship(
        "Folder", backref=backref("special_elements", cascade="all, delete-orphan")
    )

    def __repr__(self):
        return f"<SpecialElement '{self.name}' folder_id={self.folder_id}>"


class AdditionalDocumentXml(Base):
    __tablename__ = "additional_document_xml"

    document_urn = Column(
        ForeignKey("annotated_document.urn", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    )
    element_key = Column(String(255), primary_key=True, nullable=False)
    g_user_id = Column(
        ForeignKey("g_user.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    )
    xml_element = Column(String, nullable=False)
    geneious_major_version_1 = Column(
        Integer, primary_key=True, nullable=False, server_default=text("0:::INT8")
    )
    geneious_major_version_2 = Column(
        Integer, primary_key=True, nullable=False, server_default=text("0:::INT8")
    )

    g_user = relationship("GUser")

    def __repr__(self):
        return f"<AdditionalDocumentXml {self.document_urn} element_key={self.element_key}>"


class AdditionalXmlTimestamp(Base):
    __tablename__ = "additional_xml_timestamp"

    document_urn = Column(
        ForeignKey("annotated_document.urn", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    )
    g_user_id = Column(
        ForeignKey("g_user.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    )
    modified = Column(DateTime)

    g_user = relationship("GUser")

    def __repr__(self):
        return f"<AdditionalXmlTimestamp {self.document_urn} modified={self.modified}>"


class BooleanSearchFieldValue(Base):
    __tablename__ = "boolean_search_field_value"

    id = Column(
        Integer, primary_key=True, default=text(f"({NEXT_ID_SUBQUERY} {__tablename__})")
    )
    annotated_document_id = Column(
        ForeignKey("annotated_document.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )
    search_field_code = Column(
        ForeignKey("search_field.code", ondelete="CASCADE"), nullable=False, index=True
    )
    value = Column(Boolean, nullable=False)

    annotated_document = relationship(
        "AnnotatedDocument", backref=backref("boolean_search_field_values")
    )
    search_field = relationship(
        "SearchField",
        backref=backref("boolean_search_field_values", cascade="all, delete-orphan"),
    )


class DateSearchFieldValue(Base):
    __tablename__ = "date_search_field_value"

    id = Column(
        Integer, primary_key=True, default=text(f"({NEXT_ID_SUBQUERY} {__tablename__})")
    )
    annotated_document_id = Column(
        ForeignKey("annotated_document.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )
    search_field_code = Column(
        ForeignKey("search_field.code", ondelete="CASCADE"), nullable=False, index=True
    )
    value = Column(Date, nullable=False)

    annotated_document = relationship(
        "AnnotatedDocument",
        backref=backref("date_search_field_values", cascade="all, delete-orphan"),
    )
    search_field = relationship(
        "SearchField",
        backref=backref("date_search_field_values", cascade="all, delete-orphan"),
    )


t_document_read = Table(
    "document_read",
    metadata,
    Column(
        "g_user_id",
        ForeignKey("g_user.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    ),
    Column(
        "annotated_document_id",
        ForeignKey("annotated_document.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    ),
)


t_document_to_file_data = Table(
    "document_to_file_data",
    metadata,
    Column(
        "document_urn",
        ForeignKey("annotated_document.urn", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    ),
    Column(
        "file_data_id",
        ForeignKey("document_file_data.id", ondelete="CASCADE"),
        primary_key=True,
        nullable=False,
        index=True,
    ),
)


class DoubleSearchFieldValue(Base):
    __tablename__ = "double_search_field_value"

    id = Column(
        Integer, primary_key=True, default=text(f"({NEXT_ID_SUBQUERY} {__tablename__})")
    )
    annotated_document_id = Column(
        ForeignKey("annotated_document.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )
    search_field_code = Column(
        ForeignKey("search_field.code", ondelete="CASCADE"), nullable=False, index=True
    )
    value = Column(Float, nullable=False)

    annotated_document = relationship(
        "AnnotatedDocument",
        backref=backref("double_search_field_values", cascade="all, delete-orphan"),
    )
    search_field = relationship(
        "SearchField",
        backref=backref("double_search_field_values", cascade="all, delete-orphan"),
    )


class FloatSearchFieldValue(Base):
    __tablename__ = "float_search_field_value"

    id = Column(
        Integer, primary_key=True, default=text(f"({NEXT_ID_SUBQUERY} {__tablename__})")
    )
    annotated_document_id = Column(
        ForeignKey("annotated_document.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )
    search_field_code = Column(
        ForeignKey("search_field.code", ondelete="CASCADE"), nullable=False, index=True
    )
    value = Column(Float, nullable=False)

    annotated_document = relationship(
        "AnnotatedDocument",
        backref=backref("float_search_field_values", cascade="all, delete-orphan"),
    )
    search_field = relationship(
        "SearchField",
        backref=backref("float_search_field_values", cascade="all, delete-orphan"),
    )


class IntegerSearchFieldValue(Base):
    __tablename__ = "integer_search_field_value"

    id = Column(Integer, primary_key=True)
    annotated_document_id = Column(
        ForeignKey("annotated_document.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )
    search_field_code = Column(
        ForeignKey("search_field.code", ondelete="CASCADE"), nullable=False, index=True
    )
    value = Column(Integer, nullable=False)

    annotated_document = relationship(
        "AnnotatedDocument",
        backref=backref("integer_search_field_values", cascade="all, delete-orphan"),
    )
    search_field = relationship(
        "SearchField",
        backref=backref("integer_search_field_values", cascade="all, delete-orphan"),
    )


class LongSearchFieldValue(Base):
    __tablename__ = "long_search_field_value"

    id = Column(
        Integer, primary_key=True, default=text(f"({NEXT_ID_SUBQUERY} {__tablename__})")
    )
    annotated_document_id = Column(
        ForeignKey("annotated_document.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )
    search_field_code = Column(
        ForeignKey("search_field.code", ondelete="CASCADE"), nullable=False, index=True
    )
    value = Column(Integer, nullable=False)

    annotated_document = relationship(
        "AnnotatedDocument",
        backref=backref("long_search_field_values", cascade="all, delete-orphan"),
    )
    search_field = relationship(
        "SearchField",
        backref=backref("long_search_field_values", cascade="all, delete-orphan"),
    )


class StringSearchFieldValue(Base):
    __tablename__ = "string_search_field_value"

    id = Column(
        Integer, primary_key=True, default=text(f"({NEXT_ID_SUBQUERY} {__tablename__})")
    )
    annotated_document_id = Column(
        ForeignKey("annotated_document.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )
    search_field_code = Column(
        ForeignKey("search_field.code", ondelete="CASCADE"), nullable=False, index=True
    )
    value = Column(String, nullable=False)

    annotated_document = relationship(
        "AnnotatedDocument",
        backref=backref("string_field_search_values", cascade="all, delete-orphan"),
    )
    search_field = relationship(
        "SearchField",
        backref=backref("string_field_search_values", cascade="all, delete-orphan"),
    )
